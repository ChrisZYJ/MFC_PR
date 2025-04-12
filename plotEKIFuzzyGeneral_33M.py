import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.animation import FuncAnimation
import matplotlib as mpl
import re

# ==============================================
# INPUT PARAMETERS - Edit these as needed
# ==============================================

# Log file path
log_file_path = 'log/test33M.log'

# Extract case label from log file path
case_label = re.search(r'test(\w+)\.log', log_file_path).group(1)

# Define truth values
truth = {
    'G': 1e9,
    'x_center': 0.5,
    'y_center': 0.5,
    'radius': 0.25,
    'a2': 0.08,
    'a3': 0.05,
    'a4': -0.04,
    'a5': -0.03,
    'a6': 0.02,
    'a30': 0.02  # Add high frequency mode to truth
}

# Define global parameter list (change this list only to update the parameters everywhere)
parameters = ['G', 'x_center', 'y_center', 'radius', 'a2', 'a3', 'a4', 'a5', 'a6']

# Define helper lists for use in drawing the shape
translation_params = ['x_center', 'y_center', 'radius']
perturb_params = sorted([p for p in parameters if p.startswith('a')], key=lambda p: int(p[1:]))

# Updated initial values
initial_means = {
    'G': 1.1e9,
    'x_center': 0.52,
    'y_center': 0.52,
    'radius': 0.24,
    'a2': 0.0,
    'a3': 0.0,
    'a4': 0.0,
    'a5': 0.0,
    'a6': 0.0
}

initial_spreads = {
    'G': 2.0e8,
    'x_center': 0.05,
    'y_center': 0.05,
    'radius': 0.02,
    'a2': 0.1,
    'a3': 0.07,
    'a4': 0.06,
    'a5': 0.04,
    'a6': 0.03
}

# ==============================================
# CODE IMPLEMENTATION - Don't edit below unless needed
# ==============================================

# Function to parse log file and extract parameters
def parse_log_file(log_file_path):
    """Parse the log file and extract the mean and spread values for each iteration."""
    with open(log_file_path, 'r') as f:
        content = f.read()
    
    # Regex pattern to capture the iteration number, spread list, and mean list
    iteration_pattern = r'Iteration (\d+), spread = \[(.*?)\], mean = \[(.*?)\]'
    matches = re.findall(iteration_pattern, content, re.DOTALL)
    
    # Sort matches by iteration number
    matches.sort(key=lambda match: int(match[0]))
    
    max_iteration = int(matches[-1][0]) if matches else 0
    
    # Initialize dictionaries for means and spreads (each is a numpy array)
    param_means = {p: np.zeros(max_iteration + 1) for p in parameters}
    param_spreads = {p: np.zeros(max_iteration + 1) for p in parameters}
    
    # Set iteration 0 values from our initial guesses
    for p in parameters:
        param_means[p][0] = initial_means[p]
        param_spreads[p][0] = initial_spreads[p]
    
    # Fill in the parsed data for iterations 1+
    for match in matches:
        iteration = int(match[0])
        idx = iteration  # iteration 0 is already set
        spread_str = re.sub(r'\s+', ' ', match[1]).strip()
        mean_str   = re.sub(r'\s+', ' ', match[2]).strip()
        
        # Convert strings to float values (expecting one value per parameter)
        spread_values = [float(x) for x in spread_str.split()]
        mean_values = [float(x) for x in mean_str.split()]
        
        if len(mean_values) != len(parameters) or len(spread_values) != len(parameters):
            raise ValueError(f"Iteration {iteration}: Expected {len(parameters)} values but got "
                             f"{len(mean_values)} means and {len(spread_values)} spreads")
        
        for i, p in enumerate(parameters):
            param_means[p][idx] = mean_values[i]
            param_spreads[p][idx] = spread_values[i]
    
    # Print the initial and first parsed values for verification
    print("Iteration 0 (Initial) values:")
    for p in parameters:
        print(f"  {p}: mean = {param_means[p][0]}, spread = {param_spreads[p][0]}")
    
    print("\nFirst parsed iteration values:")
    for p in parameters:
        print(f"  {p}: mean = {param_means[p][1]}, spread = {param_spreads[p][1]}")
    
    total_iterations = max_iteration + 1
    print(f"\nTotal iterations: {total_iterations}")
    
    return total_iterations, param_means, param_spreads

# Parse the log file
iterations, param_means, param_spreads = parse_log_file(log_file_path)

# Create a grid for plotting the fuzzy shape
n_points = 200
x = np.linspace(0, 1, n_points)
y = np.linspace(0, 1, n_points)
X, Y = np.meshgrid(x, y)

# Create a custom colormap (from transparent white to blue)
colors = [(1, 1, 1, 0), (0, 0, 1, 1)]
cmap_name = 'fuzzy_blue'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=100)

# Set up the figure and subplots for the animation
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
fig.subplots_adjust(wspace=0.3)

# Create secondary y-axis for G
ax2_twin = ax2.twinx()

# Precompute the true boundary using truth parameters
theta_truth = np.linspace(0, 2*np.pi, 1000)
perturb_truth = 0
for p in perturb_params:
    perturb_truth += truth[p] * np.sin(int(p[1:]) * theta_truth)

# Add the high frequency mode (a30) to the true boundary only
perturb_truth += truth['a30'] * np.sin(30 * theta_truth)

r_truth = truth['radius'] + perturb_truth
x_boundary_truth = truth['x_center'] + r_truth * np.cos(theta_truth)
y_boundary_truth = truth['y_center'] + r_truth * np.sin(theta_truth)

# Function to update the animation at each iteration
def update(frame):
    # Clear all axes and any previous annotations
    fig.texts.clear()
    ax1.clear()
    ax2.clear()
    ax2_twin.clear()
    
    print(f"Processing frame {frame + 1} of {iterations}")
    
    # Sample parameters for the current iteration (for drawing and for Monte Carlo)
    current_params = {p: param_means[p][frame] for p in parameters}
    current_spreads = {p: param_spreads[p][frame] for p in parameters}
    
    # Monte Carlo sampling: use translation parameters for the base position and size
    base_x_center = current_params['x_center']
    base_y_center = current_params['y_center']
    base_radius = max(0.01, current_params['radius'])
    
    n_samples = 200  # Fixed sample size per frame
    inside_count = np.zeros_like(X)
    
    for i in range(n_samples):
        # Re-sample parameters for this Monte Carlo sample
        sample_params = {p: np.random.normal(param_means[p][frame], param_spreads[p][frame]) for p in parameters}
        sample_x_center = sample_params['x_center']
        sample_y_center = sample_params['y_center']
        sample_radius = max(0.01, sample_params['radius'])
        
        # Compute the distance and angle from each grid point to the sampled center
        dx = X - sample_x_center
        dy = Y - sample_y_center
        distance = np.sqrt(dx**2 + dy**2)
        phi = np.arctan2(dy, dx)
        
        # Accumulate contributions from all perturbation parameters (e.g. a2, a3, ...)
        perturb = 0
        for p in perturb_params:
            freq = int(p[1:])
            perturb += sample_params[p] * np.sin(freq * phi)
        radius_boundary = sample_radius + perturb
        
        inside_count[distance <= radius_boundary] += 1
    
    # Convert counts to a probability map
    probability = inside_count / n_samples
    
    # Draw the fuzzy shape using imshow
    ax1.imshow(probability, extent=[0, 1, 0, 1], origin='lower', cmap=cm, vmin=0, vmax=1)
    
    # Draw the "mean" boundary using mean values for the perturb parameters
    theta = np.linspace(0, 2*np.pi, 1000)
    perturb_mean = 0
    for p in perturb_params:
        freq = int(p[1:])
        perturb_mean += param_means[p][frame] * np.sin(freq * theta)
    r_mean = param_means['radius'][frame] + perturb_mean
    x_boundary = param_means['x_center'][frame] + r_mean * np.cos(theta)
    y_boundary = param_means['y_center'][frame] + r_mean * np.sin(theta)
    ax1.plot(x_boundary, y_boundary, 'r--', linewidth=2, label='Mean Boundary')
    
    # Draw the true boundary for comparison
    ax1.plot(x_boundary_truth, y_boundary_truth, 'g-', linewidth=1.5, label='True Boundary')
    
    if frame == 0:
        ax1.set_title('Initial Guess: Shape Estimation')
    else:
        ax1.set_title(f'Iteration {frame}: Shape Estimation')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.grid(alpha=0.3)
    ax1.legend(loc='upper right')
    
    # Plot parameters as absolute values
    geom_params = translation_params + perturb_params
    iterations_array = np.arange(0, frame + 1)
    
    # Plot all parameters except G on left y-axis
    for p in geom_params:
        color = plt.cm.tab10(geom_params.index(p) % 10)  # Get a color for this parameter
        
        # Plot actual values
        ax2.plot(iterations_array, param_means[p][:frame+1], label=p, color=color)
        
        # Plot error bars
        ax2.errorbar(iterations_array, param_means[p][:frame+1], 
                    yerr=param_spreads[p][:frame+1], fmt='o',
                    markersize=4, alpha=0.7, capsize=5, elinewidth=2, color=color)
        
        # Add a dashed line for the true value
        ax2.axhline(truth[p], color=color, linestyle='--', alpha=0.6)
    
    # Plot normalized G on the right y-axis
    g_color = 'black'  # A distinct color for G
    norm_g_values = param_means['G'][:frame+1] / truth['G']
    norm_g_spreads = param_spreads['G'][:frame+1] / truth['G']
    
    # Plot G values
    g_line, = ax2_twin.plot(iterations_array, norm_g_values, label='G (normalized)', color=g_color, linestyle='-')
    
    # Plot G error bars
    ax2_twin.errorbar(iterations_array, norm_g_values, 
                     yerr=norm_g_spreads, fmt='o',
                     markersize=4, alpha=0.7, capsize=5, elinewidth=2, color=g_color)
    
    # Add a dashed line for the true G value (normalized = 1.0)
    ax2_twin.axhline(1.0, color=g_color, linestyle='--', alpha=0.6)
    
    # Set labels and limits
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Parameter Values')
    ax2_twin.set_ylabel('Normalized G Value (G/G_true)', color=g_color)
    ax2_twin.tick_params(axis='y', labelcolor=g_color)
    
    ax2.set_xlim(0, iterations-1)
    
    # Set reasonable y-limits for the geometric parameters
    y_min = min(0, min([param_means[p][:frame+1].min() - param_spreads[p][:frame+1].max() for p in geom_params]))
    y_max = max(0.6, max([param_means[p][:frame+1].max() + param_spreads[p][:frame+1].max() for p in geom_params]))
    ax2.set_ylim(y_min, y_max)
    
    # Set reasonable y-limits for normalized G
    g_min = max(0, (norm_g_values.min() - norm_g_spreads.max()))
    g_max = norm_g_values.max() + norm_g_spreads.max()
    ax2_twin.set_ylim(min(g_min, 0.8), max(g_max, 1.3))
    
    # Add legends - we need to handle both axes
    lines, labels = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper right')
    
    ax2.set_title('Parameter Convergence')
    ax2.grid(alpha=0.3)
    
    # Display relative errors as text
    error_text_parts = []
    for p in parameters:
        rel_error = abs(param_means[p][frame] / truth[p] - 1) * 100
        error_text_parts.append(f"{p}={rel_error:.2f}%")
    error_text = "Relative Errors: " + ", ".join(error_text_parts)
    plt.figtext(0.5, 0.05, error_text, ha="center", fontsize=10,
                bbox={"facecolor": "lightgray", "alpha": 1, "pad": 5})
    
    # With blit=False, we don't need to return a list of artists
    return

# Create the animation
frame_durations = [1600]            # Longer duration for iteration 0
frame_durations.extend([800] * (iterations - 1))
ani = FuncAnimation(fig, update, frames=iterations, blit=False, interval=800)

plt.tight_layout(rect=[0, 0.1, 1, 1])  # Adjust to leave room for the text at the bottom

# Output filename with case label
output_filename = f'enki_shape_evolution_{case_label}.gif'

# Save the animation using PillowWriter
from matplotlib.animation import PillowWriter
writer = PillowWriter(fps=1.25)
ani.save(output_filename, writer=writer, dpi=120)

print(f"Animation saved as '{output_filename}'")