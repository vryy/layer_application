from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *

def ExtractMaterialInfo(prop):
    mat_info = {}
    mat_info['id'] = prop.Id
    mat_info['name'] = "Material " + str(prop.Id)
    mat_info['data'] = {}

    if prop[CONSTITUTIVE_LAW].__class__.__name__ == "Isotropic3D":
        mat_info['type'] = "isotropic elastic"
        mat_info['data']['E'] = prop[YOUNG_MODULUS]
        mat_info['data']['v'] = prop[POISSON_RATIO]

    return mat_info

def ExtractElementsByPropertiesId(elements, prop_id):
    list_elem = []
    for elem in elements:
        if elem.Properties.Id == prop_id:
            list_elem.append(elem)
    return list_elem

def ExtractConditionsByName(conditions, cond_name):
    list_cond = []
    for cond in conditions:
        if cond.__class__.__name__ == cond_name:
            list_cond.append(cond)
    return list_cond
