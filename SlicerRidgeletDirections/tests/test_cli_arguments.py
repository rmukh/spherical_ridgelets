import importlib.util
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


EXTENSION_DIR = Path(__file__).resolve().parents[1]
CLI_DIR = EXTENSION_DIR / "RidgeletDirections" / "CLI"
CLI_SCRIPT = CLI_DIR / "RidgeletDirectionsCompute.py"
CLI_XML = CLI_DIR / "RidgeletDirectionsCompute.xml"

spec = importlib.util.spec_from_file_location("ridgelet_directions_cli", CLI_SCRIPT)
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)


class RidgeletDirectionsCliTests(unittest.TestCase):
    def test_cli_descriptor_is_valid_xml(self):
        ET.parse(CLI_XML)

    def test_all_descriptor_parameters_are_accepted(self):
        descriptor = ET.parse(CLI_XML).getroot()
        parser_options = {
            option
            for action in cli.build_argument_parser()._actions
            for option in action.option_strings
        }

        for parameter in descriptor.findall("./parameters/*"):
            name = parameter.findtext("name")
            if name:
                with self.subTest(parameter=name):
                    self.assertIn(f"--{name}", parser_options)

    def test_slicer_parameter_names_map_to_native_arguments(self):
        args = cli.build_argument_parser().parse_args(
            [
                "--input_dmri", "input.nrrd",
                "--sph_J", "3",
                "--output_fiber_max_ridgelets", "directions.nrrd",
                "--fista_tolerance", "0.002",
                "--is_compress",
            ]
        )

        params = {
            cli.PARAM_INPUT_DMRI: args.input_dmri,
            cli.PARAM_SPH_J: args.sph_J,
            cli.PARAM_OUTPUT_FIBER_MAX_RIDGLET: args.output_fiber_max_ridgelets,
            cli.PARAM_FISTA_TOLERANCE: args.fista_tolerance,
            cli.PARAM_IS_COMPRESS: args.is_compress,
        }

        self.assertEqual(
            cli.build_cli_arguments(params),
            ["-i", "input.nrrd", "-sj", "3", "-ft", "0.002", "-c", "-omd_r", "directions.nrrd"],
        )

    def test_standalone_option_spellings_remain_supported(self):
        args = cli.build_argument_parser().parse_args(
            ["--input-dmri", "input.nrrd", "--output-fiber-max-ridgelets", "directions.nrrd"]
        )

        self.assertEqual(args.input_dmri, "input.nrrd")
        self.assertEqual(args.output_fiber_max_ridgelets, "directions.nrrd")


if __name__ == "__main__":
    unittest.main()