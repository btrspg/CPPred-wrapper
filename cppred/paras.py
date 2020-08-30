#! /usr/bin/env python3


import pkg_resources


def get_model_range_hexamer(species='Human'):
    '''

    '''
    if species == 'Human':
        prefix = 'Human'
    elif species == 'Integrated':
        prefix = 'Integrated'
    else:
        raise ValueError('species can only be human or integrated, not ' + species)

    m_f = pkg_resources.resource_dir(__name__, 'data/{prefix}.model'.format(prefix=prefix))
    r_f = pkg_resources.resource_dir(__name__, 'data/{prefix}.range'.format(prefix=prefix))
    h_f = pkg_resources.resource_dir(__name__, 'data/{prefix}_Hexamer.tsv'.format(prefix=prefix))
    return m_f, r_f, h_f
