
def load_sevenn(config):
    from febench.calculator.sevenn_calculator import return_calc
    calc = return_calc(config)
    return calc

def load_mace(config):
    from febench.calculator.mace_calculator import return_calc
    calc = return_calc(config)
    return calc

def load_orb(config):
    from febench.calculator.orb_calculator import return_calc
    calc = return_calc(config)
    return calc

def load_dpa(config):
    from febench.calculator.dpa_calculator import return_calc
    calc = return_calc(config)
    return calc

def load_uma(config):
    from febench.calculator.uma_calculator import return_calc
    calc = return_calc(config)
    return calc

def load_esen(config):
    from febench.calculator.esen_calculator import return_calc
    calc = return_calc(config)
    return calc


def load_test(config):
    from febench.calculator.test_calculator import return_calc
    calc = return_calc(config)
    return calc



def load_calc(config):
    calc_type = config['calculator']['calc']
    if calc_type in ['omni', 'ompa']:
        calc = load_sevenn(config)

    elif calc_type == 'mace':
        calc = load_mace(config)
    
    elif calc_type == 'uma':
        calc = load_uma(config)

    elif calc_type == 'esen':
        calc = load_esen(config)

    elif calc_type == 'test':
        calc = load_test(config)

    elif calc_type == 'orb':
        calc = load_orb(config)

    elif calc_type == 'dpa':
        calc = load_dpa(config)

    return calc
