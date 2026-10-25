# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

# a dictionary with all the variable abbreviations in the appropriate categories for plev/sfc, six_hourly, 24 hour instanteous.

s2s_variables = {
        "pressure": {
            "instantaneous_parameters": ["gh","t","u","v","q","w"]
            },
        "potential" : {"potential_temp_level": ["pv"]},
        "single_level": {
            "instantaneous_6hrly":["10u","10v","mx2t6","mn2t6",],
            "averaged_24hrs":["cape","skt","sd","rsn","asn","sm20","sm100","st20","st100","2t","2d","wtmp","ci","tcc","tcw",],
            "accumulated_24hrs":["sf","ttr","slhf","ssr","str","sshf","ssrd","strd","cp","nsss","ewss","ro","sro"],
            "instantaneous_24hrs":["sp","msl","lsm","orog","slt","TC_TRACKS","MJO"],
            "accumulated_6hrly":["tp"]
            }
        }

webAPI_params={
        "gh":"156",
        "t":"130",
        "u":"131",
        "v":"132",
        "q":"133",
        "w":"135",
        "pv":"60",
        "10u":"165",
        "10v":"166",
        "mx2t6":"121",
        "mn2t6":"122",
        "cape":"59",
        "skt":"235",
        "sd":"228141",
        "rsn":"33",
        "asn":"228032",
        "sm20":"228086",
        "sm100":"228087",
        "st20":"228095",
        "st100":"228096",
        "2t":"167",
        "2d":"168",
        "wtmp":"34",
        "ci":"31",
        "tcc":"228164",
        "tcw":"136",
        "sf":"228144",
        "ttr":"179",
        "slhf":"147",
        "ssr":"176",
        "str":"177",
        "sshf":"146",
        "ssrd":"169",
        "strd":"175",
        "cp":"228143",
        "nsss":"181",
        "ewss":"180",
        "ro":"228205",
        "sro":"174008",
        "sp":"134",
        "msl":"151",
        "lsm":"172",
        "orog":"228002",
        "slt":"43",
        "tp":"228228",
        "TC_TRACKS":"N/A"
        }

ECDS_varnames={}
