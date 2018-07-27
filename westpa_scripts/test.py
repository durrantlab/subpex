with open('pcoord.txt', 'r') as file:
    lines = file.readlines()

with open('pcoord.txt', 'w') as file:
    line = lines[0][:5] + " " + lines[1][:5]
    file.write(line)