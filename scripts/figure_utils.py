
def get_pretty_species_name(species_name, include_number=False):
        
    items = species_name.split("_")
    
    pretty_name = "%s %s" % (items[0], items[1])
    
    if pretty_name == "Faecalibacterium prausnitzii":
        include_number = True
        
    if include_number:
        pretty_name += (" (%s)" % (items[2]))
   
    ## amended nomenclature to reflect Garcia-Lopez 2019
    if pretty_name == "Bacteroides vulgatus":
        pretty_name = "Phocaeicola vulgatus"
    if pretty_name == "Bacteroides massiliensis":
        pretty_name = "Phocaeicola massiliensis"
    if pretty_name == "Eubacterium eligens":
        pretty_name = "Lachnospira eligens"
        
    return pretty_name
    
def get_abbreviated_species_name(species_name):
    
    items = species_name.split("_")
    
    pretty_name = "%s. %s" % (items[0][0], items[1])

    ## amended nomenclature to reflect Garcia-Lopez 2019
    if pretty_name == "B. vulgatus":
        pretty_name = "P. vulgatus"
    if pretty_name == "B. massiliensis":
        pretty_name = "P. massiliensis"
    if pretty_name == "E. eligens":
        pretty_name = "L. eligens"    
    
    return pretty_name