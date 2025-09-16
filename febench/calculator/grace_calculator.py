from tensorpotential.calculator import grace_fm

CALC_DCT = {
    "MP_GRACE_1L_r6_4Nov2024": "GRACE-1L-MP-r6",
    "MP_GRACE_2L_r5_4Nov2024": "GRACE-2L-MP-r5",
    "MP_GRACE_2L_r6_11Nov2024": "GRACE-2L-MP-r6",
    "GRACE-1L-OAM_2Feb25": "GRACE-1L-OAM",
    "GRACE_2L_OAM_28Jan25": "GRACE-2L-OAM",
    # shortname alisases: short name to full name
    "GRACE-2L-OMAT-L": "GRACE-2L-OMAT-large-ft-E",
    "GRACE-2L-OAM-L": "GRACE-2L-OMAT-large-ft-AM",
    "GRACE-2L-OMAT-M": "GRACE-2L-OMAT-medium-ft-E",
    "GRACE-2L-OAM-M": "GRACE-2L-OMAT-medium-ft-AM",
    "GRACE-1L-OMAT-L": "GRACE-1L-OMAT-large-ft-E",
    "GRACE-1L-OAM-L": "GRACE-1L-OMAT-large-ft-AM",
    "GRACE-1L-OMAT-M": "GRACE-1L-OMAT-medium-ft-E",
    "GRACE-1L-OAM-M": "GRACE-1L-OMAT-medium-ft-AM",
}

def return_calc(config):
    conf = config['calculator']
    model = conf['model'] # grace
    modal = conf['modal'] # oam/omat

    if modal == 'oam':
        model_name = "GRACE-2L-OAM-L"
    else:
        model_name = "GRACE-2L-OMAT-L"


    model_path =  CALC_DCT[model_name]
    print(f"[GRACE] model alias={model_name}, modal(task_name)={modal}")
    print(f"[GRACE] full model name: {model_path}")


    calc = grace_fm(model_path)

    return calc
