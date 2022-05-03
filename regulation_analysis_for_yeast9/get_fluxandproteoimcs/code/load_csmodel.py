def load_csmodel():
    import cobra
    modelCN115 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelCN115.xml')
    modelCN30 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelCN30.xml')
    modelCN50 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelCN50.xml')
    modelgln_Nlim_01 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelgln_Nlim_01.xml')
    modelile_Nlim_01 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelile_Nlim_01.xml')
    modelNH4_Nlim_005 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelNH4_Nlim_005.xml')
    modelNH4_Nlim_01 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelNH4_Nlim_01.xml')
    modelNH4_Nlim_013 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelNH4_Nlim_013.xml')
    modelNH4_Nlim_018 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelNH4_Nlim_018.xml')
    modelNH4_Nlim_030 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelNH4_Nlim_030.xml')
    modelNH4_Nlim_035 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelNH4_Nlim_035.xml')
    modelphe_Nlim_01 = cobra.io.read_sbml_model(
        r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\output\modelphe_Nlim_01.xml')
    return modelCN115, modelCN50, modelCN30, modelphe_Nlim_01, \
           modelile_Nlim_01, modelNH4_Nlim_035, modelNH4_Nlim_030, \
           modelNH4_Nlim_018, modelNH4_Nlim_013, modelNH4_Nlim_01, \
           modelNH4_Nlim_005, modelgln_Nlim_01
