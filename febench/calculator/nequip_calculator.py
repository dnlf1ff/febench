from nequip.ase import NequIPCalculator

CALC_DCT = {
        'oam': '/home/jinvk/mir-group__NequIP-OAM-L__0.1.nequip.pth'
        }


def return_calc(config):
    conf = config['calculator']
    model = conf['model']
    modal = conf['modal']
    model_name = 'NequIP-OAM-L'
    model_path = CALC_DCT[modal]

    print(f"[NequIP] model={model_name}, modal(task_name)={modal}")
    print(f"[NequIP] potential path: {model_path}")

    calc = NequIPCalculator.from_compiled_model(compile_path=model_path,
                                                device='cuda')

    return calc


