import os
import sys
import json

def get_project_root() -> str:
    """
    Returns the project root by going up one level from the script's location,
    or from the current working directory if running in a notebook.
    """
    if hasattr(sys, '_getframe'):
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
        except NameError:
            script_dir = os.getcwd()  # fallback for notebooks
    else:
        script_dir = os.getcwd()

    # Go up TWO levels to get project root
    return os.path.dirname(os.path.dirname(script_dir))
    return os.path.dirname(script_dir)  # go up one level to project root


class Path:
    @staticmethod
    def exports(name: str) -> str:
        """
        Returns the full path to a file inside the 'exports' folder,
        relative to the project root.
        """
        return os.path.join(get_project_root(), 'exports', name)


class Logger:
    @staticmethod
    def cprint(message: str, color: str = '', bold: bool = False) -> None:
        """Log a message to the terminal with optional color and bold"""
        style = color
        if bold:
            style += Colors.BOLD
        print(f"{style}{message}{Colors.ENDC}")

    @staticmethod
    def header(message: str, color: str = '') -> None:
        Logger.cprint("-"*100, color)
        Logger.cprint(message, color, True)
        Logger.cprint("-"*100, color)
        
        
        #print(f"Writing to {out_path}")
        #os.makedirs(os.path.dirname(out_path), exist_ok=True)
        #with open(out_path+"curves.json", "w") as f:
        #    json.dump(CFq.toJSON(), f, separators=(",", ":"))
        #except Exception as e:
        #    Logger.cprint(f"Failed for F_{q}: {e}", Colors.FAIL)

class Data:
    @staticmethod
    def saveJSON(path: str, fileName:str, data: dict, readable=True) -> None:
        """Save a dictionary as a JSON file."""
        out_path = f"{path}/"
        try:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with open(out_path + fileName, "w") as f:

                if readable:
                    json_str = json.dumps(data, indent=2)
                    import re
                    json_str = re.sub(r'\[\s+(\d+)\s+\]', r'[\1]', json_str)
                    # Keep simple two-element arrays on one line
                    json_str = re.sub(r'\[\s+(\d+),\s+(\d+)\s+\]', r'[\1, \2]', json_str)
                    f.write(json_str)
                # Format JSON more compactly while keeping it readable
                # json_str = json.dumps(data, indent=2)
                else:

                    json.dump(data, f, separators=(",", ":"))
                # Keep simple single-element arrays on one line
                # import re
                # json_str = re.sub(r'\[\s+(\d+)\s+\]', r'[\1]', json_str)
                # Keep simple two-element arrays on one line
                # json_str = re.sub(r'\[\s+(\d+),\s+(\d+)\s+\]', r'[\1, \2]', json_str)
                # f.write(json_str)
        except Exception as e:
            Logger.cprint(f"Failed to save JSON to {path}: {e}", Colors.FAIL)

    @staticmethod
    def loadJSON(path: str) -> dict:
        """Load a dictionary from a JSON file."""
        try:
            with open(path, "r") as f:
                return json.load(f)
        except Exception as e:
            Logger.cprint(f"Failed to load JSON from {path}: {e}", Colors.FAIL)
            return {}
        Logger.cprint("-"*100, color)


class Colors:
    # Bright colors
    HEADER = "\033[95m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_WHITE = "\033[97m"

    # Standard colors
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    TURQUOISE = "\033[36m"  # Cyan/Turquoise
    WHITE = "\033[37m"

    # Named aliases
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    SUCCESS = "\033[92m"
    INFO = "\033[94m"

    # Styles
    BOLD = "\033[1m"
    DIM = "\033[2m"
    ITALIC = "\033[3m"
    UNDERLINE = "\033[4m"

    # Reset
    ENDC = "\033[0m"


class Config:
    """Global runtime configuration, set once from CLI args."""
    rank_method: str = "mod_poly"  # "auto" | "div_poly" | "mod_poly" | "invariants"
    use_true_height: bool = False
