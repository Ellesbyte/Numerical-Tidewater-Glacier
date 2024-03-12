# Ice Flow Simulation
This Python script simulates the flow of ice over a specified domain using numerical methods and physical principles. The simulation tracks the evolution of ice thickness and surface elevation over time.

# Features
- Physical Constants: The script incorporates physical constants such as gravity, ice density, and temperature to model the flow of ice accurately.
- Flow Law: The flow law governs the rate at which ice flows under different conditions, accounting for factors such as temperature and density.
- Bottom Topography: The bottom topography of the ice sheet is modeled using parameters such as slope and Gaussian distribution to create a more realistic representation.
- Surface Mass Balance: The script accounts for the net accumulation or ablation of snow and ice at the surface, influencing the evolution of the ice sheet.
- Water-Ice Interface: Conditions are implemented to detect the water-ice interface, ensuring proper handling of situations where ice thickness reaches water depth.
- Visualization: The simulation provides real-time visualization using Matplotlib, allowing users to observe the evolution of the ice sheet over time.

# Usage
1) Dependencies: Ensure you have Python installed on your system along with the necessary libraries: NumPy and Matplotlib.
2) Running the Script: Execute the Python script ice_flow_simulation.py in your preferred Python environment.
3) Observing the Simulation: The script will generate real-time visualizations of the ice sheet evolution. You can observe the changes in ice thickness and surface elevation over the course of the simulation.
4) Interpreting the Output: Pay attention to the plotted graphs, which depict the ice sheet's behavior over time. Analyze the trends and patterns to understand the dynamics of ice flow under different conditions.

# Customization
- You can adjust parameters such as grid spacing, time steps, and physical constants within the script to explore different scenarios and model configurations.
- Experiment with different bottom topography profiles and surface mass balance patterns to simulate various ice sheet environments.
