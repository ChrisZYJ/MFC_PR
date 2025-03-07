import os
import re
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures
import logging

# ----------------------------
# CONFIGURATION
# ----------------------------

# Logging configuration
logger = logging.getLogger("EKI")
logger.setLevel(logging.DEBUG)
log_dir = os.path.join(os.path.expanduser("~/source/MFC_EKI/EKI_2D/test1"))
os.makedirs(log_dir, exist_ok=True)
fh = logging.FileHandler(os.path.join(log_dir, "EKI.log"))
fh.setLevel(logging.DEBUG)
fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(ch)

# Executable and directory paths
PRE_PROCESS_EXEC = os.path.expanduser("~/source/MFC_EKI/build/install/83d76986fd/bin/pre_process")
SIMULATION_EXEC = os.path.expanduser("~/source/MFC_EKI/build/install/877bac7bd7/bin/simulation")
BASE_DIR = os.path.expanduser("~/source/MFC_EKI/EKI_2D/test1")
NUM_INSTANCES = 8

# Simulation settings
SIM_CONFIG = {
    "dt": 5e-06,
    "t_step_start": 0,
    "t_step_stop": 1000,
    "sample_rate": 10
}
t_full = np.arange(SIM_CONFIG["t_step_start"],
                   (SIM_CONFIG["t_step_stop"] + 1) * SIM_CONFIG["dt"],
                   SIM_CONFIG["dt"])

# Ensemble Kalman Inversion (EKI) settings
EKI_CONFIG = {
    "N_ensemble": 50,
    "max_iter": 20,
    "tol": 1e-4,
    "alpha": 0.4,
    "cov_inflation": 1.2  # set to 1.0 for no inflation
}
noise_level = 0.01

# New: Number of probes to use in the simulation output
NUM_PROBES = 9

# ----------------------------
# Parameter Configuration
# ----------------------------
param_names = ["G", "x_center", "x_length", "y_center", "y_length"]
PARAM_CONFIG = {
    "G": {
        "true": 1e9,
        "guess_mean": 1.5e9,
        "guess_std": 5e8,
        "clip": (1e8, np.inf),
        "update_clip": (1e8, None),
        "files": {
            "pre_process.inp": "fluid_pp(2)%G",
            "simulation.inp": "fluid_pp(2)%G"
        }
    },
    "x_center": {
        "true": 0.5,
        "guess_mean": 0.6,
        "guess_std": 0.1,
        "clip": (0.2, 0.8),
        "update_clip": (0.2, 0.8),
        "files": {
            "pre_process.inp": "patch_icpp(2)%x_centroid"
        }
    },
    "x_length": {
        "true": 0.5,
        "guess_mean": 0.4,
        "guess_std": 0.1,
        "clip": (0.1, 0.8),
        "update_clip": (0.1, 0.8),
        "files": {
            "pre_process.inp": "patch_icpp(2)%length_x"
        }
    },
    "y_center": {
        "true": 0.5,
        "guess_mean": 0.6,
        "guess_std": 0.1,
        "clip": (0.2, 0.8),
        "update_clip": (0.2, 0.8),
        "files": {
            "pre_process.inp": "patch_icpp(2)%y_centroid"
        }
    },
    "y_length": {
        "true": 0.5,
        "guess_mean": 0.4,
        "guess_std": 0.1,
        "clip": (0.1, 0.8),
        "update_clip": (0.1, 0.8),
        "files": {
            "pre_process.inp": "patch_icpp(2)%length_y"
        }
    }
}

# ----------------------------
# HELPER FUNCTIONS
# ----------------------------

def run_simulation(sim_dir):
    """Run pre_process and simulation executables, then load the probe files."""
    D_path = os.path.join(sim_dir, "D")
    if os.path.isdir(D_path):
        try:
            shutil.rmtree(D_path)
            logger.debug(f"Deleted 'D' folder in {sim_dir}")
        except Exception as e:
            logger.error(f"Failed to delete 'D' folder in {sim_dir}: {e}")
            raise e
    try:
        logger.debug(f"Running pre_process in {sim_dir}")
        subprocess.run([PRE_PROCESS_EXEC], cwd=sim_dir, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        logger.debug(f"Running simulation in {sim_dir}")
        subprocess.run([SIMULATION_EXEC], cwd=sim_dir, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logger.error(f"Simulation failed in {sim_dir}: {e}")
        raise RuntimeError(f"Simulation failed in {sim_dir}: {e}")
    
    # Loop over the expected number of probes and load each file
    probe_data_list = []
    for i in range(1, NUM_PROBES + 1):
        probe_file = os.path.join(sim_dir, "D", f"probe{i}_prim.dat")
        if not os.path.exists(probe_file):
            msg = f"Output probe file not found: {probe_file}"
            logger.error(msg)
            raise FileNotFoundError(msg)
        p_data = np.loadtxt(probe_file)
        logger.debug(f"Read probe file {probe_file} with shape {p_data.shape}")
        probe_data_list.append(p_data)
    
    # Concatenate probe data along the time axis (flattening into one long vector)
    p_full = np.concatenate(probe_data_list)
    logger.debug(f"Concatenated probe data shape: {p_full.shape}")
    return p_full

def modify_input_files(sim_dir, params):
    """
    Update input files with current parameter values.
    For each file (pre_process.inp and simulation.inp), only parameters with a defined file pattern are updated.
    """
    logger.debug(f"Modifying input files in {sim_dir} with parameters: {params}")
    for inp_file in ["pre_process.inp", "simulation.inp"]:
        file_path = os.path.join(sim_dir, inp_file)
        with open(file_path, 'r') as f:
            content = f.read()
        for key, cfg in PARAM_CONFIG.items():
            files_map = cfg.get("files", {})
            if inp_file in files_map:
                pattern_str = files_map[inp_file]
                fmt = cfg.get("fmt", ".6e")
                pattern = rf"({re.escape(pattern_str)}\s*=\s*)\S+"
                content = re.sub(pattern,
                                 lambda m: m.group(1) + f"{params[key]:{fmt}}",
                                 content)
        with open(file_path, 'w') as f:
            f.write(content)
        logger.debug(f"Modified {inp_file} in {sim_dir}")

def prepare_simulation_folder(base_folder, folder_name):
    """Create a simulation folder by copying input files from the 'truth' folder."""
    new_sim_dir = os.path.join(base_folder, folder_name)
    os.makedirs(new_sim_dir, exist_ok=True)
    truth_dir = os.path.join(BASE_DIR, "truth")
    for file_name in ["pre_process.inp", "simulation.inp"]:
        shutil.copy(os.path.join(truth_dir, file_name), os.path.join(new_sim_dir, file_name))
    logger.debug(f"Prepared simulation folder {new_sim_dir}")
    return new_sim_dir

def simulate_member(member_id, param_vector, iter_folder, expected_length):
    """Run simulation for one ensemble member and return the subsampled output."""
    sim_dir = prepare_simulation_folder(iter_folder, f"member_{member_id+1}")
    params = {key: param_vector[i] for i, key in enumerate(param_names)}
    modify_input_files(sim_dir, params)
    try:
        p_model_full = run_simulation(sim_dir)
    except Exception as e:
        logger.error(f"Simulation error for member {member_id+1}: {e}")
        p_model_full = np.full((expected_length,), np.nan)
    # Subsample the concatenated probe data
    p_model = p_model_full[::SIM_CONFIG["sample_rate"]]
    if len(p_model) > expected_length:
        p_model = p_model[:expected_length]
    elif len(p_model) < expected_length:
        p_model = np.pad(p_model, (0, expected_length - len(p_model)), mode='edge')
    logger.debug(f"Member {member_id+1} complete with output shape {p_model.shape}")
    return p_model

def enforce_patch_constraints(ensemble):
    """Ensure patch [center - length/2, center + length/2] stays within [0.1, 0.9] for both x and y dimensions."""
    for center_name, length_name in [('x_center', 'x_length'), ('y_center', 'y_length')]:
        center_idx = param_names.index(center_name)
        length_idx = param_names.index(length_name)
        for j in range(ensemble.shape[0]):
            center_val = ensemble[j, center_idx]
            length_val = ensemble[j, length_idx]
            if center_val - length_val / 2 < 0.1:
                length_val = 2 * (center_val - 0.1)
            if center_val + length_val / 2 > 0.9:
                length_val = 2 * (0.9 - center_val)
            ensemble[j, length_idx] = length_val

# ----------------------------
# MAIN SIMULATION
# ----------------------------

# Truth run: generate "observed" data with true parameters
truth_params = {key: PARAM_CONFIG[key]["true"] for key in param_names}
truth_run_dir = prepare_simulation_folder(BASE_DIR, "truth_run")
modify_input_files(truth_run_dir, truth_params)
p_full_truth = run_simulation(truth_run_dir)
# Add noise to the truth simulation and subsample the concatenated probe data
p_obs_full = p_full_truth + noise_level * np.random.randn(len(p_full_truth))
p_obs = p_obs_full[::SIM_CONFIG["sample_rate"]]
expected_length = len(p_obs)
logger.info(f"Truth run complete. Expected observation length: {expected_length}")

# Initialize ensemble using the configured guesses
N_ensemble = EKI_CONFIG["N_ensemble"]
ensemble = np.zeros((N_ensemble, len(param_names)))
for i, key in enumerate(param_names):
    cfg = PARAM_CONFIG[key]
    if key == "G":
        guesses = np.abs(np.random.normal(loc=cfg["guess_mean"], scale=cfg["guess_std"], size=N_ensemble))
    else:
        guesses = np.random.normal(loc=cfg["guess_mean"], scale=cfg["guess_std"], size=N_ensemble)
    lower, upper = cfg["clip"]
    guesses = np.clip(guesses, lower, upper)
    ensemble[:, i] = guesses
enforce_patch_constraints(ensemble)

Gamma = (noise_level ** 2) * np.eye(expected_length)
logger.info("Starting Ensemble Kalman Inversion...")

# Store convergence history
iters = []
history = {key: {"means": [], "stds": []} for key in param_names}

# Iterative EKI update with parallel execution
for it in range(EKI_CONFIG["max_iter"]):
    iter_folder = os.path.join(BASE_DIR, f"iter_{it+1}")
    os.makedirs(iter_folder, exist_ok=True)
    logger.debug(f"Iteration {it+1} in folder {iter_folder}")

    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_INSTANCES) as executor:
        results = list(executor.map(lambda j: simulate_member(j, ensemble[j, :],
                                                              iter_folder,
                                                              expected_length),
                                      range(N_ensemble)))
    Y = np.array(results)  # shape: (N_ensemble, expected_length)

    m_bar = np.mean(ensemble, axis=0)
    y_bar = np.mean(Y, axis=0)
    C_md = sum(np.outer(ensemble[j, :] - m_bar, Y[j, :] - y_bar) for j in range(N_ensemble)) / (N_ensemble - 1)
    C_dd = sum(np.outer(Y[j, :] - y_bar, Y[j, :] - y_bar) for j in range(N_ensemble)) / (N_ensemble - 1)
    K = C_md @ np.linalg.inv(C_dd + Gamma)

    # Update each ensemble member and enforce parameter constraints
    for j in range(N_ensemble):
        innovation = p_obs - Y[j, :]
        update = EKI_CONFIG["alpha"] * (K @ innovation)
        ensemble[j, :] += update
        for i, key in enumerate(param_names):
            lower, upper = PARAM_CONFIG[key]["update_clip"]
            if lower is not None and upper is not None:
                ensemble[j, i] = np.clip(ensemble[j, i], lower, upper)
            elif lower is not None:
                ensemble[j, i] = max(lower, ensemble[j, i])
            elif upper is not None:
                ensemble[j, i] = min(upper, ensemble[j, i])
    enforce_patch_constraints(ensemble)

    # Covariance inflation
    if EKI_CONFIG["cov_inflation"] != 1.0:
        m_post = np.mean(ensemble, axis=0)
        ensemble = m_post + EKI_CONFIG["cov_inflation"] * (ensemble - m_post)

    spread = np.std(ensemble, axis=0)
    mean_val = np.mean(ensemble, axis=0)
    rel_err = {key: abs(mean_val[i] - PARAM_CONFIG[key]["true"]) / PARAM_CONFIG[key]["true"] * 100
               for i, key in enumerate(param_names)}
    logger.info(f"Iteration {it+1}, spread = {spread}, mean = {mean_val}")
    logger.info("Relative errors (%): " + ", ".join(f"{key} = {rel_err[key]:.2f}%" for key in param_names))

    iters.append(it+1)
    for i, key in enumerate(param_names):
        history[key]["means"].append(mean_val[i] / PARAM_CONFIG[key]["true"])
        history[key]["stds"].append(spread[i] / PARAM_CONFIG[key]["true"])

    if np.all(spread / mean_val < EKI_CONFIG["tol"]):
        logger.info(f"Converged at iteration {it+1}")
        break

estimated_params = np.mean(ensemble, axis=0)
logger.info("Estimated parameters:")
for i, key in enumerate(param_names):
    logger.info(f"{key}: {estimated_params[i]:.6e}")

# ----------------------------
# PLOTTING RESULTS
# ----------------------------

fig, axs = plt.subplots(1, len(param_names), figsize=(20, 4))
for i, key in enumerate(param_names):
    axs[i].hist(ensemble[:, i], bins=15, alpha=0.7, label='Ensemble')
    axs[i].axvline(PARAM_CONFIG[key]["true"], color='r', linestyle='--', label='True')
    axs[i].set_xlabel(key)
    axs[i].set_ylabel('Frequency')
    axs[i].set_title(f"Distribution of {key}")
    axs[i].legend()
plt.tight_layout()
plot_filename = os.path.join(BASE_DIR, "ensemble_distribution.png")
plt.savefig(plot_filename)
plt.close()
logger.info(f"Saved ensemble distribution plots to {plot_filename}")

plt.figure(figsize=(10, 6))
markers = ['-o', '-s', '-^', '-d', '-v']
for i, key in enumerate(param_names):
    plt.errorbar(iters, history[key]["means"], yerr=history[key]["stds"],
                 fmt=markers[i], capsize=5, label=f"{key} (normalized)")
plt.xlabel("Iteration")
plt.ylabel("Normalized Parameter Value")
plt.title("Convergence History")
plt.legend()
plt.grid(True)
conv_plot_filename = os.path.join(BASE_DIR, "convergence_history.png")
plt.savefig(conv_plot_filename)
plt.close()
logger.info(f"Saved convergence history plot to {conv_plot_filename}")
