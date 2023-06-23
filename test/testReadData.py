import sys
sys.path.append("..")
import offset_simulation as off

pp, rst, rstRaw = off.read_pressure_data("../data/pressure_scan/20230519/pressure_scan_data.root", fmt='numpy')

print(pp)
print(type(pp))
