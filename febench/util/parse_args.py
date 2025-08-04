import argparse

def parse_base_args(argv: list[str]| None=None): #TODO namespace
    # -----------BASE------------#
    parser = argparse.ArgumentParser(description= "cli tool")

    parser.add_argument('--calc', type=str, default='omni',help='potential name')
    
    parser.add_argument('--modal', type=str, default='mpa',help='data modal; mpa, omat24')
    
    parser.add_argument('--config', type=str, default='./config.yaml', help='config yaml file dir')

    parser.add_argument('--carbon_config', type=str, default='./carbon_config.yaml', help='carbon config yaml file dir')

    return parser.parse_args(argv)

