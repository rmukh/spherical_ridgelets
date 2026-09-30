#!/usr/bin/env python
"""CLI module for computing spherical ridgelets direct fiber directions."""

import os
import subprocess
import sys
from typing import Optional

# CLI parameter names matching the XML
PARAM_INPUT_DMRI = "input_dmri"
PARAM_INPUT_MASK = "input_mask"
PARAM_EXTERNAL_GRADIENTS = "external_gradients"
PARAM_SPH_J = "sph_J"
PARAM_SPH_RHO = "sph_rho"
PARAM_LVL = "lvl"
PARAM_N_SPLITS = "n_splits"
PARAM_FISTA_LAMBDA = "fista_lambda"
PARAM_FISTA_ITERATIONS = "fista_iterations"
PARAM_FISTA_TOLERANCE = "fista_tolerance"
PARAM_NTH = "nth"
PARAM_IS_COMPRESS = "is_compress"
PARAM_OUTPUT_FIBER_MAX_RIDGLET = "output_fiber_max_ridgelets"
PARAM_OUTPUT_RIDGELETS = "output_ridgelets"
PARAM_OUTPUT_SIGNAL_RECON = "output_signal_recon"
PARAM_OUTPUT_ODF = "output_odf"
PARAM_RIDGLET_NMS_ANGLE = "ridgelet_nms_angle"
PARAM_MAX_ODF_THRESH = "max_odf_thresh"
PARAM_PRINT_SCALE_WEIGHTS = "print_scale_weights"
PARAM_TEST_SCALE_WEIGHTS = "test_scale_weights"
PARAM_TEST_DIRECT_RIDGLET_MAXIMA = "test_direct_ridgelet_maxima"


def _get_executable_path() -> Optional[str]:
    """Find the sphridg executable."""
    # Check common locations
    possible_paths = [
        # From Slicer extension installation
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "bin", "sphridg.exe"),
        os.path.join(os.path.dirname(__file__), "sphridg.exe"),
        # Environment variable
        os.environ.get("SPHRIDG_EXECUTABLE", ""),
    ]
    
    # Also check the build directory relative to this module
    module_dir = os.path.dirname(__file__)
    for root, dirs, files in os.walk(module_dir):
        if "sphridg.exe" in files:
            return os.path.join(root, "sphridg.exe")
    
    for path in possible_paths:
        if path and os.path.isfile(path):
            return path
    
    return None


def build_cli_arguments(params: dict) -> list:
    """Build command-line arguments for sphridg.exe from parameter dict."""
    args = []
    
    # Input dMRI volume (required)
    if params.get(PARAM_INPUT_DMRI):
        args.extend(["-i", params[PARAM_INPUT_DMRI]])
    
    # Input mask (optional)
    if params.get(PARAM_INPUT_MASK):
        args.extend(["-m", params[PARAM_INPUT_MASK]])
    
    # External gradients (optional)
    if params.get(PARAM_EXTERNAL_GRADIENTS):
        args.extend(["-ext_grads", params[PARAM_EXTERNAL_GRADIENTS]])
    
    # Spherical ridgelets parameters
    if params.get(PARAM_SPH_J) is not None:
        args.extend(["-sj", str(params[PARAM_SPH_J])])
    
    if params.get(PARAM_SPH_RHO) is not None:
        args.extend(["-srho", str(params[PARAM_SPH_RHO])])
    
    if params.get(PARAM_LVL) is not None:
        args.extend(["-lvl", str(params[PARAM_LVL])])
    
    if params.get(PARAM_N_SPLITS) is not None:
        args.extend(["-nspl", str(params[PARAM_N_SPLITS])])
    
    if params.get(PARAM_FISTA_LAMBDA) is not None:
        args.extend(["-lmd", str(params[PARAM_FISTA_LAMBDA])])
    
    if params.get(PARAM_FISTA_ITERATIONS) is not None:
        args.extend(["-fi", str(params[PARAM_FISTA_ITERATIONS])])
    
    if params.get(PARAM_FISTA_TOLERANCE) is not None:
        args.extend(["-ft", str(params[PARAM_FISTA_TOLERANCE])])
    
    if params.get(PARAM_NTH) is not None:
        args.extend(["-nth", str(params[PARAM_NTH])])
    
    if params.get(PARAM_IS_COMPRESS):
        args.append("-c")
    
    # Output options
    if params.get(PARAM_OUTPUT_FIBER_MAX_RIDGLET):
        args.extend(["-omd_r", params[PARAM_OUTPUT_FIBER_MAX_RIDGLET]])
    
    if params.get(PARAM_OUTPUT_RIDGELETS):
        args.extend(["-ridg", params[PARAM_OUTPUT_RIDGELETS]])
    
    if params.get(PARAM_OUTPUT_SIGNAL_RECON):
        args.extend(["-sr", params[PARAM_OUTPUT_SIGNAL_RECON]])
    
    if params.get(PARAM_OUTPUT_ODF):
        args.extend(["-odf", params[PARAM_OUTPUT_ODF]])
    
    if params.get(PARAM_RIDGLET_NMS_ANGLE) is not None:
        args.extend(["-rth", str(params[PARAM_RIDGLET_NMS_ANGLE])])
    
    if params.get(PARAM_MAX_ODF_THRESH) is not None:
        args.extend(["-mth", str(params[PARAM_MAX_ODF_THRESH])])
    
    # Test/debug options
    if params.get(PARAM_PRINT_SCALE_WEIGHTS):
        args.append("-sw")
    
    if params.get(PARAM_TEST_SCALE_WEIGHTS):
        args.append("-test_sw")
    
    if params.get(PARAM_TEST_DIRECT_RIDGLET_MAXIMA):
        args.append("-test_omd_r")
    
    return args


def compute_ridgelet_directions(params: dict, progress_callback=None) -> tuple:
    """
    Compute ridgelet directions by calling the sphridg executable.
    
    Returns:
        tuple: (success: bool, error_message: str or None)
    """
    exe_path = _get_executable_path()
    if not exe_path:
        return False, "Could not find sphridg.exe executable. Please set SPHRIDG_EXECUTABLE environment variable or ensure the executable is in the expected location."
    
    if not os.path.isfile(exe_path):
        return False, f"sphridg.exe not found at: {exe_path}"
    
    args = build_cli_arguments(params)
    
    # Check that at least one output is specified
    output_params = [
        PARAM_OUTPUT_FIBER_MAX_RIDGLET,
        PARAM_OUTPUT_RIDGELETS,
        PARAM_OUTPUT_SIGNAL_RECON,
        PARAM_OUTPUT_ODF,
    ]
    has_output = any(params.get(p) for p in output_params)
    has_test = params.get(PARAM_PRINT_SCALE_WEIGHTS) or params.get(PARAM_TEST_SCALE_WEIGHTS) or params.get(PARAM_TEST_DIRECT_RIDGLET_MAXIMA)
    
    if not has_output and not has_test:
        return False, "At least one output must be specified."
    
    # Build full command
    cmd = [exe_path] + args
    
    try:
        # Run the executable
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout
        )
        
        if result.returncode != 0:
            error_msg = result.stderr or result.stdout or "Unknown error"
            return False, f"sphridg.exe failed with return code {result.returncode}: {error_msg}"
        
        return True, None
        
    except subprocess.TimeoutExpired:
        return False, "Computation timed out after 1 hour."
    except Exception as e:
        return False, f"Error running sphridg.exe: {str(e)}"


def build_argument_parser():
    """Build an argument parser for Slicer and standalone CLI invocations."""
    import argparse

    parser = argparse.ArgumentParser(description="Compute spherical ridgelets direct fiber directions")

    # Input parameters
    parser.add_argument("--input-dmri", "--input_dmri", required=True, help="Input dMRI volume file")
    parser.add_argument("--input-mask", "--input_mask", help="Optional brain mask volume file")
    parser.add_argument("--external-gradients", "--external_gradients", help="Optional external gradient directions file")

    # Spherical ridgelets parameters
    parser.add_argument("--sph-J", "--sph_J", type=int, default=2, help="Spherical ridgelets J (default: 2)")
    parser.add_argument("--sph-rho", "--sph_rho", type=float, default=3.125, help="Spherical ridgelets rho (default: 3.125)")
    parser.add_argument("--lvl", type=int, default=4, help="Icosahedron tesselation order (default: 4)")
    parser.add_argument("--n-splits", "--n_splits", type=int, default=-1, help="Number of splits (default: -1 for auto)")
    parser.add_argument("--fista-lambda", "--fista_lambda", type=float, default=0.01, help="FISTA lambda (default: 0.01)")
    parser.add_argument("--fista-iterations", "--fista_iterations", type=int, default=100, help="FISTA iterations (default: 100)")
    parser.add_argument("--fista-tolerance", "--fista_tolerance", type=float, default=0.0001, help="FISTA tolerance (default: 0.0001)")
    parser.add_argument("--nth", type=int, default=-1, help="Number of threads (default: -1 for auto)")
    parser.add_argument("--is-compress", "--is_compress", action="store_true", help="Enable compression")

    # Output parameters
    parser.add_argument("--output-fiber-max-ridgelets", "--output_fiber_max_ridgelets", help="Output direct ridgelet maxima file (-omd_r)")
    parser.add_argument("--output-ridgelets", "--output_ridgelets", help="Output ridgelets coefficients file")
    parser.add_argument("--output-signal-recon", "--output_signal_recon", help="Output signal reconstruction file")
    parser.add_argument("--output-odf", "--output_odf", help="Output ODF values file")

    # Additional parameters
    parser.add_argument("--ridgelet-nms-angle", "--ridgelet_nms_angle", type=float, default=20.0, help="Direct ridgelet NMS angle (default: 20)")
    parser.add_argument("--max-odf-thresh", "--max_odf_thresh", type=float, default=0.7, help="ODF maxima threshold (default: 0.7)")

    # Test options
    parser.add_argument("--print-scale-weights", "--print_scale_weights", action="store_true", help="Print scale weights")
    parser.add_argument("--test-scale-weights", "--test_scale_weights", action="store_true", help="Test scale weights")
    parser.add_argument("--test-direct-ridgelet-maxima", "--test_direct_ridgelet_maxima", action="store_true", help="Test direct ridgelet maxima")
    return parser


def main():
    """Main entry point for the CLI module."""
    args = build_argument_parser().parse_args()

    # Build params dict
    params = {
        PARAM_INPUT_DMRI: args.input_dmri,
        PARAM_INPUT_MASK: args.input_mask,
        PARAM_EXTERNAL_GRADIENTS: args.external_gradients,
        PARAM_SPH_J: args.sph_J,
        PARAM_SPH_RHO: args.sph_rho,
        PARAM_LVL: args.lvl,
        PARAM_N_SPLITS: args.n_splits,
        PARAM_FISTA_LAMBDA: args.fista_lambda,
        PARAM_FISTA_ITERATIONS: args.fista_iterations,
        PARAM_FISTA_TOLERANCE: args.fista_tolerance,
        PARAM_NTH: args.nth,
        PARAM_IS_COMPRESS: args.is_compress,
        PARAM_OUTPUT_FIBER_MAX_RIDGLET: args.output_fiber_max_ridgelets,
        PARAM_OUTPUT_RIDGELETS: args.output_ridgelets,
        PARAM_OUTPUT_SIGNAL_RECON: args.output_signal_recon,
        PARAM_OUTPUT_ODF: args.output_odf,
        PARAM_RIDGLET_NMS_ANGLE: args.ridgelet_nms_angle,
        PARAM_MAX_ODF_THRESH: args.max_odf_thresh,
        PARAM_PRINT_SCALE_WEIGHTS: args.print_scale_weights,
        PARAM_TEST_SCALE_WEIGHTS: args.test_scale_weights,
        PARAM_TEST_DIRECT_RIDGLET_MAXIMA: args.test_direct_ridgelet_maxima,
    }
    
    success, error = compute_ridgelet_directions(params)
    
    if success:
        print("Ridgelet directions computed successfully.")
        sys.exit(0)
    else:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

