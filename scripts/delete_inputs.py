import sys
from pathlib import Path

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python bash_test.py <var1>")
        sys.exit(1)

    filename = "inputs/inputs_" + sys.argv[1] # string
    
    if not filename.endswith(".txt"):
        filename = filename + ".txt"

    file_path = Path(filename)

    if file_path.exists():
        file_path.unlink()
    else:
        print(f"{file_path} does not exist")
