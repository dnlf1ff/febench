
def load_sevenn(config):
    from febench.calculator.sevenn_calculator import return_calc
    calc = return_calc(config)
    return calc

def load_mace(config):
    from febench.calculator.mace_calculator import return_calc
    calc = return_calc(config)
    return calc

def load_uma(config):
    from febench.calculator.uma_calculator import return_calc
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

    elif calc_type == 'test':
        calc = load_test(config)

    return calc
