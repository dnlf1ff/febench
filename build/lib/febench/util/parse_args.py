import argparse

def parse_base_args(argv: list[str]| None=None): #TODO namespace
    parser = argparse.ArgumentParser(description= "cli tool")

    parser.add_argument('--config', type=str, default='./config.yaml', 
                        help='config yaml file directory')

    parser.add_argument('--calc_type', type=str, default='sevennet-mf',
                        help='sevennet,sevennet-mf, esen, uma, etc')
 
    parser.add_argument('--prefix', '--cwd', dest='prefix', type=str, default='omni-omat24',
                        help='working directory')
    
    parser.add_argument('--model','--calc', dest='model', type=str, default='.',
                        help='absolute directory of the U-MLIP potential')

    parser.add_argument('--model_dirname','--model_path','--path','--dirname', dest='dirname',type=str, default='.',
                        help='absolute directory of the U-MLIP potential')


    parser.add_argument('--modal','--task', '--head','--task_name', dest='modal', type=str, 
                        default='omat24',help='data modality for multi-functional U-MLIPs; mpa, omat24 etc.')

    parser.add_argument('--dispersion', type=bool, 
                        default=False, help='whether to exclude D3')

    parser.add_argument('--functional', type=str, 
                        default='PBE', help='Functional')


    return parser.parse_args(argv)

