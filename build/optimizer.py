#!/usr/bin/env python3
"""
Complete optimization example for your muon tracking setup
CSV format: Detector,ParticleType,x[cm],px[MeV/c],y[cm],py[MeV/c],z[cm],pz[MeV/c],TotalEnergy[MeV]
"""

import os
import sys
import numpy as np
from solenoid_optimizer import *
from muon_csv_analyzer import MuonTrackingAnalyzer

class MuonLossOptimizer:
    """Simplified optimizer specifically for your muon tracking setup"""
    
    def __init__(self, executable_path: str, gap_bounds: Tuple[float, float] = (5.0, 100.0)):
        self.executable_path = executable_path
        self.gap_bounds = gap_bounds
        self.results_history = []
        
        # Create workspace
        self.workspace = "muon_optimization"
        os.makedirs(self.workspace, exist_ok=True)
    
    def run_single_simulation(self, gap1: float, gap2: float, n_runs: int = 1000) -> Dict:
        """Run a single simulation with specified gaps"""
        
        # Create unique run directory
        run_id = len(self.results_history) + 1
        run_dir = os.path.join(self.workspace, f"run_{run_id:03d}")
        os.makedirs(run_dir, exist_ok=True)
        
        # Output CSV file
        output_csv = os.path.join(run_dir, "particles.csv")
        
        # Run your executable
        # Assuming command format: ./executable gap1 gap2 n_runs output_file
        cmd = f"{self.executable_path} {gap1} {gap2} {n_runs} {output_csv}"
        
        print(f"Running: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=run_dir)
        
        if result.returncode != 0:
            raise RuntimeError(f"Simulation failed: {result.stderr}")
        
        # Analyze results
        if not os.path.exists(output_csv):
            raise FileNotFoundError(f"Expected output file not found: {output_csv}")
        
        # Parse CSV using analyzer
        analyzer = MuonTrackingAnalyzer(output_csv)
        metrics = analyzer.calculate_optimization_metrics()
        
        # Store results
        simulation_result = {
            'run_id': run_id,
            'gaps': [gap1, gap2],
            'n_runs': n_runs,
            'detector_counts': metrics['detector_counts'],
            'transmission_rates': metrics['transmission_rates'],
            'total_loss': metrics['total_loss'],
            'beam_quality': metrics['beam_quality'],
            'output_file': output_csv
        }
        
        self.results_history.append(simulation_result)
        
        return simulation_result
    
    def objective_function(self, gaps: List[float]) -> float:
        """
        Objective function for optimization
        Returns a score to minimize (lower is better)
        """
        try:
            result = self.run_single_simulation(gaps[0], gaps[1])
            
            # Multi-objective scoring
            total_loss_weight = 0.7
            beam_quality_weight = 0.3
            
            # Primary objective: minimize total particle loss
            loss_score = result['total_loss']
            
            # Secondary objective: beam quality (normalize spreads)
            quality_score = 0
            if result['beam_quality']:
                # Normalize quality metrics (adjust scaling based on your typical values)
                energy_spread = result['beam_quality'].get('energy_spread', 0) / 10.0  # Assume ~10 MeV typical
                spatial_spread = (
                    result['beam_quality'].get('spatial_spread_x', 0) + 
                    result['beam_quality'].get('spatial_spread_y', 0)
                ) / 2.0 / 5.0  # Assume ~5 cm typical
                
                quality_score = (energy_spread + spatial_spread) / 2.0
            
            # Combined score
            combined_score = total_loss_weight * loss_score + beam_quality_weight * quality_score
            
            print(f"Gaps: [{gaps[0]:.1f}, {gaps[1]:.1f}] -> Loss: {loss_score:.3f}, Quality: {quality_score:.3f}, Score: {combined_score:.3f}")
            
            return combined_score
            
        except Exception as e:
            print(f"Simulation failed for gaps {gaps}: {e}")
            return 1.0  # Return worst possible score
    
    def optimize_bayesian(self, n_evaluations: int = 30) -> Dict:
        """Run Bayesian optimization"""
        from skopt import gp_minimize
        from skopt.space import Real
        
        print(f"Starting Bayesian optimization with {n_evaluations} evaluations")
        print(f"Gap bounds: {self.gap_bounds}")
        
        # Define search space
        space = [
            Real(self.gap_bounds[0], self.gap_bounds[1], name='gap1'),
            Real(self.gap_bounds[0], self.gap_bounds[1], name='gap2')
        ]
        
        # Run optimization
        result = gp_minimize(
            func=self.objective_function,
            dimensions=space,
            n_calls=n_evaluations,
            n_initial_points=8,
            acq_func='EI',  # Expected Improvement
            random_state=42
        )
        
        # Extract best result
        best_gaps = result.x
        best_score = result.fun
        
        # Run final verification with more statistics
        print(f"\nVerifying optimal solution...")
        final_result = self.run_single_simulation(best_gaps[0], best_gaps[1], n_runs=5000)
        
        optimization_result = {
            'best_gaps': best_gaps,
            'best_score': best_score,
            'final_verification': final_result,
            'optimization_history': self.results_history,
            'skopt_result': result
        }
        
        return optimization_result
    
    def grid_search(self, n_points_per_dim: int = 8) -> Dict:
        """Simple grid search for comparison/validation"""
        gap1_values = np.linspace(self.gap_bounds[0], self.gap_bounds[1], n_points_per_dim)
        gap2_values = np.linspace(self.gap_bounds[0], self.gap_bounds[1], n_points_per_dim)
        
        print(f"Running grid search: {n_points_per_dim}x{n_points_per_dim} = {n_points_per_dim**2} evaluations")
        
        best_score = float('inf')
        best_gaps = None
        grid_results = []
        
        total_evaluations = n_points_per_dim * n_points_per_dim
        current_eval = 0
        
        for gap1 in gap1_values:
            for gap2 in gap2_values:
                current_eval += 1
                print(f"Grid search progress: {current_eval}/{total_evaluations}")
                
                score = self.objective_function([gap1, gap2])
                grid_results.append({
                    'gaps': [gap1, gap2],
                    'score': score
                })
                
                if score < best_score:
                    best_score = score
                    best_gaps = [gap1, gap2]
        
        # Final verification
        final_result = self.run_single_simulation(best_gaps[0], best_gaps[1], n_runs=5000)
        
        return {
            'best_gaps': best_gaps,
            'best_score': best_score,
            'final_verification': final_result,
            'grid_results': grid_results
        }
    
    def analyze_optimization_results(self, optimization_result: Dict):
        """Analyze and display optimization results"""
        best_gaps = optimization_result['best_gaps']
        final_result = optimization_result['final_verification']
        
        print("\n" + "="*80)
        print("OPTIMIZATION RESULTS SUMMARY")
        print("="*80)
        
        print(f"Optimal Gap Configuration:")
        print(f"  Gap 1 (Section 1→2): {best_gaps[0]:.2f} mm")
        print(f"  Gap 2 (Section 2→3): {best_gaps[1]:.2f} mm")
        
        print(f"\nDetector Performance:")
        detector_counts = final_result['detector_counts']
        transmission_rates = final_result['transmission_rates']
        
        for i, (det_id, count) in enumerate(detector_counts.items()):
            rate = transmission_rates[i] if i < len(transmission_rates) else 0
            print(f"  Detector {det_id}: {count} particles ({rate:.3f} transmission)")
        
        print(f"\nOverall Performance:")
        print(f"  Total Particle Loss: {final_result['total_loss']:.3f}")
        print(f"  Final Transmission: {transmission_rates[-1]:.3f}" if transmission_rates else "  No transmission data")
        
        if final_result['beam_quality']:
            print(f"\nBeam Quality Metrics:")
            for metric, value in final_result['beam_quality'].items():
                print(f"  {metric}: {value:.3f}")
        
        # Calculate improvement over worst case
        if len(self.results_history) > 1:
            all_losses = [r['total_loss'] for r in self.results_history]
            worst_loss = max(all_losses)
            improvement = (worst_loss - final_result['total_loss']) / worst_loss * 100
            print(f"\nImprovement over worst configuration: {improvement:.1f}% reduction in loss")
    
    def plot_optimization_progress(self, optimization_result: Dict):
        """Plot optimization progress and results"""
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Muon Loss Optimization Results')
        
        # Optimization progress
        losses = [r['total_loss'] for r in self.results_history]
        axes[0,0].plot(losses, 'b-o', markersize=4)
        axes[0,0].set_xlabel('Evaluation Number')
        axes[0,0].set_ylabel('Total Loss')
        axes[0,0].set_title('Optimization Progress')
        axes[0,0].grid(True)
        
        # Gap parameter space (if we have enough points)
        if len(self.results_history) >= 10:
            gap1_vals = [r['gaps'][0] for r in self.results_history]
            gap2_vals = [r['gaps'][1] for r in self.results_history]
            colors = [r['total_loss'] for r in self.results_history]
            
            scatter = axes[0,1].scatter(gap1_vals, gap2_vals, c=colors, cmap='viridis_r', s=50)
            axes[0,1].set_xlabel('Gap 1 [mm]')
            axes[0,1].set_ylabel('Gap 2 [mm]')
            axes[0,1].set_title('Parameter Space Exploration')
            plt.colorbar(scatter, ax=axes[0,1], label='Total Loss')
            
            # Mark best point
            best_gaps = optimization_result['best_gaps']
            axes[0,1].scatter(best_gaps[0], best_gaps[1], c='red', s=200, marker='*', 
                            edgecolor='black', linewidth=2, label='Optimal')
            axes[0,1].legend()
        
        # Transmission rates comparison
        if 'final_verification' in optimization_result:
            final_result = optimization_result['final_verification']
            transmission_rates = final_result['transmission_rates']
            detectors = list(range(1, len(transmission_rates) + 1))
            
            axes[1,0].bar(detectors, transmission_rates, alpha=0.7)
            axes[1,0].set_xlabel('Detector Number')
            axes[1,0].set_ylabel('Transmission Rate')
            axes[1,0].set_title('Final Transmission Rates')
            axes[1,0].set_ylim(0, 1)
            
            # Add loss annotations
            for i, rate in enumerate(transmission_rates):
                if i > 0:
                    loss = transmission_rates[i-1] - rate
                    axes[1,0].text(i+1, rate + 0.02, f'Loss: {loss:.3f}', 
                                  ha='center', fontsize=8)
        
        # Beam quality metrics (if available)
        if 'final_verification' in optimization_result and optimization_result['final_verification']['beam_quality']:
            quality = optimization_result['final_verification']['beam_quality']
            metrics = list(quality.keys())
            values = list(quality.values())
            
            axes[1,1].bar(range(len(metrics)), values, alpha=0.7)
            axes[1,1].set_xticks(range(len(metrics)))
            axes[1,1].set_xticklabels([m.replace('_', '\n') for m in metrics], rotation=45, ha='right')
            axes[1,1].set_ylabel('Value')
            axes[1,1].set_title('Beam Quality Metrics')
        else:
            axes[1,1].text(0.5, 0.5, 'No beam quality data', ha='center', va='center', transform=axes[1,1].transAxes)
        
        plt.tight_layout()
        plt.show()
    
    def save_results(self, optimization_result: Dict, filename: str = None):
        """Save optimization results to file"""
        if filename is None:
            filename = f"muon_optimization_results.json"
        
        filepath = os.path.join(self.workspace, filename)
        
        # Convert numpy arrays and complex objects to JSON-serializable format
        import json
        
        def convert_for_json(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif hasattr(obj, '__dict__'):
                return str(obj)
            else:
                return obj
        
        # Create serializable copy
        serializable_result = {}
        for key, value in optimization_result.items():
            if key == 'skopt_result':
                # Skip the scikit-optimize result object for now
                continue
            else:
                serializable_result[key] = convert_for_json(value)
        
        with open(filepath, 'w') as f:
            json.dump(serializable_result, f, indent=2)
        
        print(f"Results saved to: {filepath}")
        
        # Also save a summary text file
        summary_file = filepath.replace('.json', '_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("MUON LOSS OPTIMIZATION SUMMARY\n")
            f.write("="*50 + "\n\n")
            f.write(f"Optimal gaps: {optimization_result['best_gaps']}\n")
            f.write(f"Best score: {optimization_result['best_score']:.6f}\n")
            if 'final_verification' in optimization_result:
                final = optimization_result['final_verification']
                f.write(f"Final loss: {final['total_loss']:.6f}\n")
                f.write(f"Transmission rates: {final['transmission_rates']}\n")

def main():
    """Main execution function - customize this for your setup"""
    
    # CONFIGURATION - MODIFY THESE PATHS FOR YOUR SETUP
    executable_path = "./your_muon_simulation_executable"  # Path to your GEANT4 executable
    gap_bounds = (10.0, 80.0)  # Min and max gap sizes in mm
    
    # Check if executable exists
    if not os.path.exists(executable_path):
        print(f"ERROR: Executable not found at {executable_path}")
        print("Please update the executable_path in the main() function")
        return
    
    # Initialize optimizer
    optimizer = MuonLossOptimizer(executable_path, gap_bounds)
    
    print("Starting muon loss optimization...")
    print(f"Executable: {executable_path}")
    print(f"Gap bounds: {gap_bounds} mm")
    
    # Choose optimization method
    method = input("Choose optimization method (b=Bayesian, g=Grid): ").lower().strip()
    
    if method.startswith('g'):
        # Grid search
        n_points = int(input("Grid points per dimension (default 6): ") or "6")
        results = optimizer.grid_search(n_points)
    else:
        # Bayesian optimization (default)
        n_evals = int(input("Number of evaluations (default 25): ") or "25")
        results = optimizer.optimize_bayesian(n_evals)
    
    # Analyze results
    optimizer.analyze_optimization_results(results)
    
    # Plot results
    try:
        optimizer.plot_optimization_progress(results)
    except ImportError:
        print("Matplotlib not available for plotting")
    
    # Save results
    optimizer.save_results(results)
    
    return results

def quick_test():
    """Quick test function to verify your setup works"""
    executable_path = "./your_muon_simulation_executable"  # Update this
    
    if not os.path.exists(executable_path):
        print(f"Please update executable_path in quick_test() function")
        return
    
    optimizer = MuonLossOptimizer(executable_path, gap_bounds=(20.0, 60.0))
    
    print("Running quick test with gaps [30, 40] mm...")
    try:
        result = optimizer.run_single_simulation(30.0, 40.0, n_runs=500)
        print(f"Test successful!")
        print(f"  Transmission rates: {result['transmission_rates']}")
        print(f"  Total loss: {result['total_loss']:.3f}")
        return True
    except Exception as e:
        print(f"Test failed: {e}")
        return False

if __name__ == "__main__":
    # Uncomment one of these options:
    
    # Option 1: Quick test to verify setup
    # quick_test()
    
    # Option 2: Run full optimization
    results = main()
    
    # Option 3: Analyze existing CSV file
    # csv_file = "path/to/your/output.csv"
    # analyzer = MuonTrackingAnalyzer(csv_file)
    # analyzer.generate_summary_report()
    # analyzer.plot_detector_comparison()

"""
USAGE INSTRUCTIONS:

1. Update the executable_path in main() function to point to your GEANT4 executable

2. Make sure your executable accepts command line arguments in this format:
   ./your_executable gap1_mm gap2_mm n_runs output_csv_file

3. Your executable should output a CSV file with this exact format:
   Detector,ParticleType,x[cm],px[MeV/c],y[cm],py[MeV/c],z[cm],pz[MeV/c],TotalEnergy[MeV]

4. Run the optimization:
   python muon_optimization_example.py

5. Choose optimization method:
   - Bayesian (b): More efficient, good for expensive simulations
   - Grid (g): Systematic exploration, good for validation

6. Results will be saved in the 'muon_optimization' directory

EXAMPLE COMMAND LINE INTERFACE FOR YOUR EXECUTABLE:
Your GEANT4 program should accept:
./your_executable 25.5 33.2 1000 /path/to/output.csv

Where:
- 25.5 = gap1 in mm
- 33.2 = gap2 in mm  
- 1000 = number of simulation runs
- /path/to/output.csv = output file path
"""