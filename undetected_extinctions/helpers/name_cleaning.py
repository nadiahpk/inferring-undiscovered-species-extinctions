import re

# a list of things that I will do to each species name to make them consistent etc.
def clean_species_name(spp_name):

    # remove unicode \xa0s, non-breaking space in Latin1 (ISO 8859-1), also chr(160), with a space
    spp_name = spp_name.replace(u'\xa0', u' ')

    # make the whole name lower case 
    spp_name = spp_name.lower()

    # then uppercase the first letter
    #spp_name = spp_name[:1].upper() + spp_name[1:]

    # remove double or more spaces
    spp_name = re.sub(' +', ' ', spp_name)

    # remove brackets, just makes my life easier
    spp_name = spp_name.replace('(','').replace(')','')

    # strip beginning and ending whitespaces
    spp_name = spp_name.strip()

    return spp_name
