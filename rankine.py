from steam import steam


class rankine:
    """
    A class to model the Rankine power cycle, which converts heat into work.
    The efficiency of the cycle is determined by the properties of steam at various points.
    """

    def __init__(self, p_low, p_high, t_high=None, name='Rankine Cycle'):
        """
        Initializes the Rankine cycle with specified parameters.

        :param p_low: The low pressure of the cycle in kPa.
        :param p_high: The high pressure of the cycle in kPa.
        :param t_high: The temperature at the high pressure, if superheated steam is used.
        :param name: A name for the cycle for identification.
        """
        self.p_low = p_low
        self.p_high = p_high
        self.t_high = t_high
        self.name = name
        self.efficiency = None
        self.turbine_work = None
        self.pump_work = None
        self.heat_added = None

        self.state1 = None
        self.state2 = None
        self.state3 = None
        self.state4 = None

    def calc_efficiency(self):
        """
        Calculates the efficiency of the Rankine cycle and the associated thermodynamic properties.
        """
        if self.t_high is None:
            self.state1 = steam(self.p_high, x=1, name='Turbine Inlet')
        else:
            self.state1 = steam(self.p_high, T=self.t_high, name='Turbine Inlet')
        self.state1.calc()

        self.state2 = steam(self.p_low, s=self.state1.s, name='Turbine Exit')
        self.state2.calc()

        self.state3 = steam(self.p_low, x=0, name='Pump Inlet')
        self.state3 = steam(self.p_low, x=0, name='Pump Inlet')  # Saturated liquid at low pressure
        self.state3.calc()

        # Check if state3's properties, particularly 'h' and 'v', were correctly set
        if self.state3.h is None or self.state3.v is None:
            print("Debug info for state 3:")
            print(f"h: {self.state3.h}, v: {self.state3.v}")
            raise ValueError("State 3 properties could not be calculated. Check input data and steam table.")

        # Assuming state3's properties are now correctly calculated, proceed with state4 calculation
        self.state4 = steam(self.p_high, s=self.state3.s, name='Pump Exit')
        self.state4.calc()
        self.state4.h = self.state3.h + self.state3.v * (self.p_high - self.p_low) / 1000

        # Continue with efficiency calculation...
        self.turbine_work = self.state1.h - self.state2.h
        self.pump_work = self.state4.h - self.state3.h
        self.heat_added = self.state1.h - self.state4.h

        self.efficiency = (self.turbine_work - self.pump_work) / self.heat_added * 100

    def print_summary(self):
        """
        Prints a summary of the Rankine cycle's efficiency and state properties.
        """
        if self.efficiency is None:
            self.calc_efficiency()

        print(f'Cycle Summary for: {self.name}')
        print(f'\tEfficiency: {self.efficiency:.2f}%')
        print(f'\tTurbine Work: {self.turbine_work:.2f} kJ/kg')
        print(f'\tPump Work: {self.pump_work:.2f} kJ/kg')
        print(f'\tHeat Added: {self.heat_added:.2f} kJ/kg\n')

        for state in [self.state1, self.state2, self.state3, self.state4]:
            state.print()


# Example usage
def main():
    cycle = rankine(p_low=8, p_high=8000, t_high=None, name='Example Rankine Cycle')
    cycle.calc_efficiency()
    cycle.print_summary()


if __name__ == "__main__":
    main()
