# test_rankine.py

from rankine import RankineCycle

# Test for saturated vapor entering turbine
saturated_cycle = RankineCycle(p_high=8000, p_low=8, x1=1)
saturated_efficiency = saturated_cycle.efficiency()
print("Efficiency of Rankine cycle with saturated vapor: {:.2f}%".format(saturated_efficiency * 100))

# Test for superheated steam entering turbine
superheated_cycle = RankineCycle(p_high=8000, p_low=8, superheat_temp=1.7 * 250)  # 250 is placeholder for Tsat at p_high
superheated_efficiency = superheated_cycle.efficiency()
print("Efficiency of Rankine cycle with superheated vapor: {:.2f}%".format(superheated_efficiency * 100))
