import json
import ConfigParser
import argparse
import datetime


def timestr():
    return datetime.datetime.strftime(datetime.datetime.now(), 
            '%Y-%m-%d-%H-%M')


def open_config(infile):
    """ creates a config parser and opens infile
    returns config object"""
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.readfp(open(infile))
    return config

def add_to_config(config, section, key, value):
    """ can add new item to a config file
    creates new section if necessary

    eg 
    section = 'subject_info'
    key = 'subid'
    value = 'despo_211'
    will create this in the config file
    [subject_info]
    subid = despo_211

    Note
    ----
    all inputs should be strings to write properly
    """
    if not config.has_section(section):
        config.add_section(section)
    config.set(section, key, value)


def config_to_dict(config):
    outd = {}
    sections = config.sections()
    for section in sections:
        options = config.options(section)
        for option in options:
            try:
                outd[option] = json.loads(config.get(section, option))
            except ValueError:
                outd[option] = config.get(section, option)
    return outd


def write_new_config(config, outfile):
    """ write new config to new outfile
    useful for updating config for subject specific values
    and writing to subject directory"""
    with open(outfile, 'wb') as configfile:
        config.write(configfile)