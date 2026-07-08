import sys

path = sys.argv[1]
with open(path, "r") as f:
    c = f.read()
c = c.replace("arg_data args[Size];", "arg_data args[Size == 0 ? 1 : Size];")
with open(path, "w") as f:
    f.write(c)
