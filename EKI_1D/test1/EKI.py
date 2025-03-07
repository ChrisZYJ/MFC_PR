import os
import re
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures
import logging

# ----------------------------
# Logging Configuration
# ----------------------------
logger = logging.getLogger("EKI")
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler(os.path.join(os.path.expanduser("~/source/MFC_EKI/EKI_1D/test1"), "EKI.log"))
fh.setLevel(logging.DEBUG)
fh_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(fh_formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch_formatter = logging.Formatter('%(message)s')
ch.setFormatter(ch_formatter)
logger.addHandler(ch)

# ----------------------------
# Fortran Executable Paths (adjust if needed)
# ----------------------------
PRE_PROCESS_EXEC = os.path.expanduser("~/source/MFC_EKI/build/install/83d76986fd/bin/pre_process")
SIMULATION_EXEC = os.path.expanduser("~/source/MFC_EKI/build/install/877bac7bd7/bin/simulation")

# ----------------------------
# Base Directory
# ----------------------------
BASE_DIR = os.path.expanduser("~/source/MFC_EKI/EKI_1D/test1")

# ----------------------------
# Parallel Execution Settings
# ----------------------------
NUM_INSTANCES = 8

# ----------------------------
# Global Simulation Parameters
# ----------------------------
dt = 5e-06
t_step_start = 0
t_step_stop = 1000
t_full = np.arange(t_step_start, (t_step_stop + 1) * dt, dt)
sample_rate = 10  # sample every 10th time step

# ----------------------------
# Helper Functions
# ----------------------------
def run_simulation(sim_dir):
    """
    Execute the Fortran forward solver in the given directory.
    Before running, delete any existing 'D' folder to avoid contaminating probe data.
    Runs pre_process then simulation and reads the probe file.
    """
    # Delete the 'D' folder if it exists (only the D folder, not the whole sim_dir)
    D_path = os.path.join(sim_dir, "D")
    if os.path.exists(D_path) and os.path.isdir(D_path):
        try:
            shutil.rmtree(D_path)
            logger.debug(f"Deleted existing 'D' folder in {sim_dir}")
        except Exception as e:
            logger.error(f"Failed to delete 'D' folder in {sim_dir}: {e}")
            raise e

    try:
        logger.debug(f"Running pre_process in {sim_dir}")
        subprocess.run([PRE_PROCESS_EXEC], cwd=sim_dir, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.debug(f"Running simulation in {sim_dir}")
        subprocess.run([SIMULATION_EXEC], cwd=sim_dir, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logger.error(f"Simulation failed in {sim_dir}: {e}")
        raise RuntimeError(f"Simulation failed in {sim_dir}: {e}")
    
    probe_file = os.path.join(sim_dir, "D", "probe1_prim.dat")
    if not os.path.exists(probe_file):
        msg = f"Output probe file not found in {sim_dir}/D/probe1_prim.dat"
        logger.error(msg)
        raise FileNotFoundError(msg)
    
    p_full = np.loadtxt(probe_file)
    logger.debug(f"Read probe file in {sim_dir} with shape {p_full.shape}")
    return p_full

def modify_input_files(sim_dir, G, center, length):
    """
    Modify the input files in sim_dir:
      - Replace fluid_pp(2)%G with G (both files).
      - For pre_process.inp, update patch_icpp(2)%x_centroid to center 
        and patch_icpp(2)%length_x to length.
    """
    logger.debug(f"Modifying input files in {sim_dir}: G={G}, center={center}, length={length}")
    
    for inp_file in ["pre_process.inp", "simulation.inp"]:
        file_path = os.path.join(sim_dir, inp_file)
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Update fluid_pp(2)%G for both files.
        content = re.sub(
            r"(fluid_pp\(2\)%G\s*=\s*)\S+",
            lambda m: m.group(1) + f"{G:.6e}",
            content
        )
        if inp_file == "pre_process.inp":
            # Set patch_icpp(2)%x_centroid = center and patch_icpp(2)%length_x = length.
            content = re.sub(
                r"(patch_icpp\(2\)%x_centroid\s*=\s*)\S+",
                lambda m: m.group(1) + f"{center:.6f}",
                content
            )
            content = re.sub(
                r"(patch_icpp\(2\)%length_x\s*=\s*)\S+",
                lambda m: m.group(1) + f"{length:.6f}",
                content
            )
        with open(file_path, 'w') as f:
            f.write(content)
        logger.debug(f"Modified {inp_file} in {sim_dir}")

def prepare_simulation_folder(base_folder, folder_name):
    """
    Create a new folder (folder_name) in base_folder by copying the input files
    from the 'truth' folder. Returns the new simulation folder path.
    """
    new_sim_dir = os.path.join(base_folder, folder_name)
    os.makedirs(new_sim_dir, exist_ok=True)
    
    truth_dir = os.path.join(BASE_DIR, "truth")
    for file_name in ["pre_process.inp", "simulation.inp"]:
        src = os.path.join(truth_dir, file_name)
        dst = os.path.join(new_sim_dir, file_name)
        shutil.copy(src, dst)
    logger.debug(f"Prepared simulation folder {new_sim_dir}")
    return new_sim_dir

def simulate_member(member_id, ensemble_params, iter_folder, expected_length):
    """
    Create simulation folder for an ensemble member, modify inputs,
    run the simulation, and return the subsampled pressure time series.
    """
    sim_folder_name = f"member_{member_id+1}"
    sim_dir = prepare_simulation_folder(iter_folder, sim_folder_name)
    G, center, length = ensemble_params
    modify_input_files(sim_dir, G, center, length)
    try:
        p_model_full = run_simulation(sim_dir)
    except Exception as e:
        logger.error(f"Simulation error for ensemble member {member_id+1}: {e}")
        p_model_full = np.full_like(t_full, np.nan)
    p_model = p_model_full[::sample_rate]
    if len(p_model) > expected_length:
        p_model = p_model[:expected_length]
    elif len(p_model) < expected_length:
        p_model = np.pad(p_model, (0, expected_length - len(p_model)), mode='edge')
    logger.debug(f"Member {member_id+1} simulation complete with output shape {p_model.shape}")
    return p_model

# ----------------------------
# Generate "Observed" Data (Truth Run)
# ----------------------------
# True parameters:
true_G = 1e9         # Truth for Young's Modulus (G)
true_center = 0.5    # Truth for patch center
true_length = 0.5    # Truth for patch length

truth_run_dir = prepare_simulation_folder(BASE_DIR, "truth_run")
modify_input_files(truth_run_dir, true_G, true_center, true_length)
p_full_truth = run_simulation(truth_run_dir)
noise_level = 0.01
p_obs_full = p_full_truth + noise_level * np.random.randn(len(p_full_truth))
p_obs = p_obs_full[::sample_rate]
expected_length = len(p_obs)
logger.info(f"Truth run complete. Expected observation length: {expected_length}")

# ----------------------------
# Ensemble Kalman Inversion Settings
# ----------------------------
N_ensemble = 50
max_iter = 20
tol = 1e-4
alpha = 0.4

ensemble = np.zeros((N_ensemble, 3))
# Generate guesses for G around 5e9 (e.g., mean 5e9, std 1e7)
ensemble[:, 0] = np.abs(np.random.normal(loc=1.5e9, scale=5e8, size=N_ensemble))
# Generate guesses for center around 0.5 (bounded between 0.4 and 0.6)
ensemble[:, 1] = np.clip(np.random.normal(loc=0.6, scale=0.15, size=N_ensemble), 0.2, 0.8)
# Generate guesses for length around 0.5 (bounded between 0.4 and 0.6)
ensemble[:, 2] = np.random.normal(loc=0.4, scale=0.15, size=N_ensemble)

# Enforce that patch remains within [0.1, 0.9]: center - length/2 >= 0.1 and center + length/2 <= 0.9.
for j in range(N_ensemble):
    center_val = ensemble[j, 1]
    length_val = ensemble[j, 2]
    if center_val - length_val/2 < 0.1:
        length_val = 2 * (center_val - 0.1)
    if center_val + length_val/2 > 0.9:
        length_val = 2 * (0.9 - center_val)
    ensemble[j, 2] = length_val
    
Gamma = (noise_level**2) * np.eye(expected_length)
logger.info("Starting Ensemble Kalman Inversion...")

# ----------------------------
# Convergence History Storage
# ----------------------------
iters = []
G_means = []
G_stds = []
center_means = []
center_stds = []
length_means = []
length_stds = []

# ----------------------------
# Iterative EKI Update with Parallel Execution
# ----------------------------
for it in range(max_iter):
    iter_folder = os.path.join(BASE_DIR, f"iter_{it+1}")
    os.makedirs(iter_folder, exist_ok=True)
    logger.debug(f"Starting iteration {it+1} in folder {iter_folder}")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_INSTANCES) as executor:
        results = executor.map(
            lambda j: simulate_member(j, ensemble[j, :], iter_folder, expected_length),
            range(N_ensemble)
        )
    Y = np.array(list(results))  # Y shape: (N_ensemble, expected_length)
    
    m_bar = np.mean(ensemble, axis=0)      # shape: (3,)
    y_bar = np.mean(Y, axis=0)             # shape: (expected_length,)
    
    C_md = np.zeros((3, expected_length))
    for j in range(N_ensemble):
        C_md += np.outer(ensemble[j, :] - m_bar, Y[j, :] - y_bar)
    C_md /= (N_ensemble - 1)
    
    C_dd = np.zeros((expected_length, expected_length))
    for j in range(N_ensemble):
        diff = Y[j, :] - y_bar
        C_dd += np.outer(diff, diff)
    C_dd /= (N_ensemble - 1)
    
    K = C_md @ np.linalg.inv(C_dd + Gamma)
    
    for j in range(N_ensemble):
        innovation = p_obs - Y[j, :]  # both shape: (expected_length,)
        update = alpha * (K @ innovation)  # shape: (3,)
        ensemble[j, :] += update
        
        # Enforce constraints:
        ensemble[j, 0] = max(1e8, ensemble[j, 0])  # G must be positive (min 1e8)
        ensemble[j, 1] = np.clip(ensemble[j, 1], 0.2, 0.8)
        ensemble[j, 2] = ensemble[j, 2]
        
        # Enforce that patch remains within [0.1, 0.9]: center - length/2 >= 0.1 and center + length/2 <= 0.9.
        for j in range(N_ensemble):
            center = ensemble[j, 1]
            length = ensemble[j, 2]
            if center - length/2 < 0.1:
                length = 2 * (center - 0.1)
            if center + length/2 > 0.9:
                length = 2 * (0.9 - center)
            ensemble[j, 2] = length
    
    spread = np.std(ensemble, axis=0)
    mean_val = np.mean(ensemble, axis=0)
    logger.info(f"Iteration {it+1}, ensemble spread = {spread}, ensemble mean = {mean_val}")
    
    # Store convergence history normalized by true parameter values.
    iters.append(it+1)
    G_means.append(mean_val[0] / true_G)
    G_stds.append(spread[0] / true_G)
    center_means.append(mean_val[1] / true_center)
    center_stds.append(spread[1] / true_center)
    length_means.append(mean_val[2] / true_length)
    length_stds.append(spread[2] / true_length)
    
    if np.all(spread/mean_val < tol):
        logger.info(f"Converged at iteration {it+1}")
        break

estimated_params = np.mean(ensemble, axis=0)
logger.info("Estimated parameters:")
logger.info(f"Young's Modulus (G): {estimated_params[0]:.6e}")
logger.info(f"Patch center: {estimated_params[1]:.6f}")
logger.info(f"Patch length: {estimated_params[2]:.6f}")

# ----------------------------
# Save the Final Ensemble Distribution Plots for each parameter
# ----------------------------
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
axs[0].hist(ensemble[:, 0], bins=15, alpha=0.7, label='Ensemble')
axs[0].axvline(true_G, color='r', linestyle='--', label='True G')
axs[0].set_xlabel("Young's Modulus (G)")
axs[0].set_ylabel('Frequency')
axs[0].set_title("Final Ensemble Distribution of G")
axs[0].legend()

axs[1].hist(ensemble[:, 1], bins=15, alpha=0.7, label='Ensemble')
axs[1].axvline(true_center, color='r', linestyle='--', label='True center')
axs[1].set_xlabel('Patch Center')
axs[1].set_ylabel('Frequency')
axs[1].set_title('Final Ensemble Distribution of Patch Center')
axs[1].legend()

axs[2].hist(ensemble[:, 2], bins=15, alpha=0.7, label='Ensemble')
axs[2].axvline(true_length, color='r', linestyle='--', label='True length')
axs[2].set_xlabel('Patch Length')
axs[2].set_ylabel('Frequency')
axs[2].set_title('Final Ensemble Distribution of Patch Length')
axs[2].legend()

plt.tight_layout()
plot_filename = os.path.join(BASE_DIR, "ensemble_distribution.png")
plt.savefig(plot_filename)
plt.close()
logger.info(f"Saved ensemble distribution plots to {plot_filename}")

# ----------------------------
# Save Convergence History Plot
# ----------------------------
plt.figure(figsize=(10, 6))
# Plot for Young's Modulus (G)
plt.errorbar(iters, G_means, yerr=G_stds, fmt='-o', label="G (normalized)")
# Plot for patch center
plt.errorbar(iters, center_means, yerr=center_stds, fmt='-s', label="Center (normalized)")
# Plot for patch length
plt.errorbar(iters, length_means, yerr=length_stds, fmt='-^', label="Length (normalized)")

plt.xlabel("Iteration")
plt.ylabel("Normalized Parameter Value")
plt.title("Convergence History of Parameters")
plt.legend()
plt.grid(True)
conv_plot_filename = os.path.join(BASE_DIR, "convergence_history.png")
plt.savefig(conv_plot_filename)
plt.close()
logger.info(f"Saved convergence history plot to {conv_plot_filename}")
