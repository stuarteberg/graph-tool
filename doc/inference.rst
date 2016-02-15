.. automodule:: graph_tool.inference
    :no-undoc-members:
    :show-inheritance:

.. testcode:: inference_detailed
   :hide:

   import test_inference

.. testoutput:: inference_detailed
   :hide:
   :options: -ELLIPSIS, +NORMALIZE_WHITESPACE

   directed: True overlap: False layered: False deg-corr: False dl: False
            mcmc (unweighted)
                    (653.7368247664258, 100) 72
            mcmc
                    (709.1284329780722, 103) 67
                    (26.555814693961594, 102) 68
                    (-108.54759577258069, 115) 73
            merge
                    (327.7692003897153, 50)
                    (240.08404807493113, 110)
                    (80.17105004298814, 115)
            shrink
                    5

   directed: True overlap: False layered: False deg-corr: False dl: True
            mcmc (unweighted)
                    (37.79906768205631, 105) 77
            mcmc
                    (55.97800577471912, 107) 72
                    (-83.84539912129875, 99) 61
                    (60.363819212130494, 114) 73
            merge
                    (-403.6228272212937, 50)
                    (-2.5333139051913243, 104)
                    (19.27758710817686, 115)
            shrink
                    5

   directed: True overlap: False layered: False deg-corr: True dl: False
            mcmc (unweighted)
                    (448.9786966541199, 102) 79
            mcmc
                    (523.3965114137179, 106) 74
                    (-25.920853037839876, 103) 74
                    (25.289255024521793, 113) 74
            merge
                    (277.5286735275392, 50)
                    (358.72935230058386, 101)
                    (-115.7860369979704, 115)
            shrink
                    5

   directed: True overlap: False layered: False deg-corr: True dl: True
            mcmc (unweighted)
                    (-370.02300429105844, 97) 68
            mcmc
                    (-236.5251583848658, 104) 77
                    (-90.66222214522463, 104) 71
                    (22.38348339446021, 112) 72
            merge
                    (-747.2417795796297, 50)
                    (-137.5051248827061, 100)
                    (-22.19343083384164, 113)
            shrink
                    5

   directed: True overlap: False layered: covariates deg-corr: False dl: False
            mcmc (unweighted)
                    (676.3594764878382, 114) 72
            mcmc
                    (659.612425515866, 109) 73
                    (-7.114556958148668, 115) 75
                    (12.450844591087305, 115) 73
            merge
                    (546.8411931837169, 50)
                    (446.7766774894063, 111)
                    (-233.62505581815967, 115)
            shrink
                    5

   directed: True overlap: False layered: covariates deg-corr: False dl: True
            mcmc (unweighted)
                    (11.134477993729128, 113) 71
            mcmc
                    (39.69116889632038, 110) 70
                    (23.949240858672248, 114) 72
                    (-20.311156765735674, 115) 75
            merge
                    (-185.8315002121304, 50)
                    (-8.659869616217907, 113)
                    (7.1698821344710115, 114)
            shrink
                    5

   directed: True overlap: False layered: covariates deg-corr: True dl: False
            mcmc (unweighted)
                    (625.1252715930277, 113) 70
            mcmc
                    (687.1284762185751, 113) 64
                    (88.17735417358503, 113) 59
                    (-267.6672513613203, 114) 73
            merge
                    (461.52365374179743, 50)
                    (554.2643481837697, 109)
                    (-311.10773578568643, 114)
            shrink
                    5

   directed: True overlap: False layered: covariates deg-corr: True dl: True
            mcmc (unweighted)
                    (-265.4943379846294, 111) 76
            mcmc
                    (-260.48845827199153, 115) 75
                    (-41.58828679977694, 114) 72
                    (20.07067071778924, 114) 74
            merge
                    (-584.6644045599711, 50)
                    (-124.7745711915025, 112)
                    (70.39505964132627, 115)
            shrink
                    5

   directed: True overlap: False layered: True deg-corr: False dl: False
            mcmc (unweighted)
                    (649.827058787757, 109) 63
            mcmc
                    (460.97415340386925, 115) 74
                    (20.20550023851644, 113) 71
                    (-11.482268963688943, 114) 75
            merge
                    (382.0710697233173, 50)
                    (659.1089055216413, 112)
                    (-384.8843394117986, 115)
            shrink
                    5

   directed: True overlap: False layered: True deg-corr: False dl: True
            mcmc (unweighted)
                    (-146.93776337069448, 114) 67
            mcmc
                    (-118.54952723454903, 112) 67
                    (-22.25121342666771, 111) 70
                    (0.09101846869118524, 113) 72
            merge
                    (-345.98550748381786, 50)
                    (-19.856093212550192, 111)
                    (24.406541323617336, 115)
            shrink
                    5

   directed: True overlap: False layered: True deg-corr: True dl: False
            mcmc (unweighted)
                    (299.47637000759045, 113) 71
            mcmc
                    (339.7882309572504, 113) 73
                    (-22.951725336949753, 113) 70
                    (-71.18214115859378, 115) 76
            merge
                    (112.65674339105865, 50)
                    (270.7282008346159, 115)
                    (-48.227031269529164, 115)
            shrink
                    5

   directed: True overlap: False layered: True deg-corr: True dl: True
            mcmc (unweighted)
                    (-489.54010405218867, 115) 72
            mcmc
                    (-477.0610489947019, 115) 73
                    (-44.315852689841975, 111) 70
                    (-6.031490138403797, 113) 70
            merge
                    (-694.3072536306113, 50)
                    (-260.61832985101154, 112)
                    (32.31080727878356, 115)
            shrink
                    5

   directed: True overlap: True layered: False deg-corr: False dl: False
            mcmc (unweighted)
                    (684.2122420163477, 1222) 784
                    (-14.493178594335856, 1340) 787
            merge
                    (0.0, 50)
            shrink
                    5

   directed: True overlap: True layered: False deg-corr: False dl: True
            mcmc (unweighted)
                    (-451.273408828854, 1223) 777
                    (-27.657518258027803, 1340) 769
            merge
                    (-308.96675893350766, 50)
            shrink
                    5

   directed: True overlap: True layered: False deg-corr: True dl: False
            mcmc (unweighted)
                    (410.11174108492406, 1226) 756
                    (-11.783502069519063, 1337) 760
            merge
                    (1.3862943611198906, 50)
            shrink
                    5

   directed: True overlap: True layered: False deg-corr: True dl: True
            mcmc (unweighted)
                    (-160.75919034773452, 1225) 785
                    (-37.94994320552587, 1340) 763
            merge
                    (-326.8859937813693, 50)
            shrink
                    5

   directed: True overlap: True layered: covariates deg-corr: False dl: False
            mcmc (unweighted)
                    (6.4407729047092115, 1197) 115
                    (15.090003687444668, 1306) 115
            merge
                    (-53.01104703350246, 50)
            shrink
                    5

   directed: True overlap: True layered: covariates deg-corr: False dl: True
            mcmc (unweighted)
                    (0.7469127360648091, 1189) 115
                    (-4.454066106756823, 1307) 115
            merge
                    (-1060.3966480302265, 50)
            shrink
                    5

   directed: True overlap: True layered: covariates deg-corr: True dl: False
            mcmc (unweighted)
                    (19.34655635911719, 1201) 115
                    (13.021627155555198, 1302) 115
            merge
                    (139.75704771951064, 50)
            shrink
                    5

   directed: True overlap: True layered: covariates deg-corr: True dl: True
            mcmc (unweighted)
                    (-42.609890899551715, 1196) 115
                    (-20.85691258865186, 1310) 115
            merge
                    (-1234.0387385376016, 50)
            shrink
                    5

   directed: True overlap: True layered: True deg-corr: False dl: False
            mcmc (unweighted)
                    (-32.53034167952721, 1201) 115
                    (56.45069492553456, 1294) 115
            merge
                    (-114.8981755030546, 50)
            shrink
                    5

   directed: True overlap: True layered: True deg-corr: False dl: True
            mcmc (unweighted)
                    (-1.035017095430904, 1191) 115
                    (-4.378896497938934, 1305) 115
            merge
                    (-993.1502579858155, 50)
            shrink
                    5

   directed: True overlap: True layered: True deg-corr: True dl: False
            mcmc (unweighted)
                    (6.508901229588357, 1189) 115
                    (-13.11412599497536, 1301) 115
            merge
                    (-77.33333924925819, 50)
            shrink
                    5

   directed: True overlap: True layered: True deg-corr: True dl: True
            mcmc (unweighted)
                    (62.56398098817939, 1195) 115
                    (-3.5377428709883647, 1304) 115
            merge
                    (-1052.7716211643888, 50)
            shrink
                    5

   directed: False overlap: False layered: False deg-corr: False dl: False
            mcmc (unweighted)
                    (693.169084122821, 106) 69
            mcmc
                    (575.1255994230743, 106) 74
                    (-18.268855666875268, 102) 74
                    (84.23562085572648, 114) 68
            merge
                    (254.062569976686, 50)
                    (294.07064838987327, 108)
                    (-35.54039844306018, 115)
            shrink
                    5

   directed: False overlap: False layered: False deg-corr: False dl: True
            mcmc (unweighted)
                    (6.089470037990972, 99) 74
            mcmc
                    (31.337138743014393, 96) 71
                    (-32.348734568160104, 94) 68
                    (11.13293108914912, 114) 73
            merge
                    (-452.35567315108005, 50)
                    (-14.94414868569813, 108)
                    (-26.47985161842624, 115)
            shrink
                    5

   directed: False overlap: False layered: False deg-corr: True dl: False
            mcmc (unweighted)
                    (474.71492598670625, 100) 79
            mcmc
                    (575.0190616712077, 103) 74
                    (35.03115802167634, 99) 72
                    (-44.92187114668639, 114) 73
            merge
                    (249.97743129324903, 50)
                    (425.23174309430004, 110)
                    (-159.32783404320577, 114)
            shrink
                    5

   directed: False overlap: False layered: False deg-corr: True dl: True
            mcmc (unweighted)
                    (-193.2001940253293, 99) 74
            mcmc
                    (-174.90080436460175, 102) 76
                    (-39.52484390424932, 104) 72
                    (48.052018152160485, 113) 72
            merge
                    (-717.2044026904009, 50)
                    (-71.77424635033832, 105)
                    (-25.114234051660894, 114)
            shrink
                    5

   directed: False overlap: False layered: covariates deg-corr: False dl: False
            mcmc (unweighted)
                    (772.5453980570721, 114) 67
            mcmc
                    (720.1691858052512, 112) 69
                    (68.31407337248119, 111) 65
                    (-146.74986338638587, 114) 72
            merge
                    (473.38852115957025, 50)
                    (376.67392866586084, 114)
                    (1.1626252512182127, 114)
            shrink
                    5

   directed: False overlap: False layered: covariates deg-corr: False dl: True
            mcmc (unweighted)
                    (39.60067070934042, 115) 72
            mcmc
                    (-9.194100392728705, 115) 77
                    (23.758390517950442, 110) 75
                    (7.399872077643064, 115) 78
            merge
                    (-254.49516803072956, 50)
                    (-9.215586113589342, 112)
                    (23.632624271547694, 112)
            shrink
                    5

   directed: False overlap: False layered: covariates deg-corr: True dl: False
            mcmc (unweighted)
                    (679.6815899057593, 111) 70
            mcmc
                    (737.5066421528536, 114) 68
                    (-41.630424206875865, 111) 68
                    (-71.64216522170558, 114) 73
            merge
                    (483.74983009336034, 50)
                    (465.972158053492, 114)
                    (-168.20407971455091, 114)
            shrink
                    5

   directed: False overlap: False layered: covariates deg-corr: True dl: True
            mcmc (unweighted)
                    (-246.9705243300452, 113) 67
            mcmc
                    (-236.51398935319358, 113) 68
                    (83.02986995537107, 113) 78
                    (-33.29065578699246, 114) 73
            merge
                    (-504.63687228651776, 50)
                    (-171.50178946195857, 108)
                    (49.50677409599346, 113)
            shrink
                    5

   directed: False overlap: False layered: True deg-corr: False dl: False
            mcmc (unweighted)
                    (540.8663037198617, 113) 68
            mcmc
                    (507.7047847062117, 115) 69
                    (-56.73592493601191, 113) 74
                    (14.478960856237514, 114) 75
            merge
                    (385.119877617005, 50)
                    (609.4426771976545, 112)
                    (-95.63740447873565, 114)
            shrink
                    5

   directed: False overlap: False layered: True deg-corr: False dl: True
            mcmc (unweighted)
                    (-135.48400339347597, 114) 69
            mcmc
                    (-177.55059658149057, 112) 67
                    (77.30739145915234, 114) 77
                    (-37.788312169996864, 115) 78
            merge
                    (-378.09608309928564, 50)
                    (-9.759576031614198, 114)
                    (-0.9536321167456725, 115)
            shrink
                    5

   directed: False overlap: False layered: True deg-corr: True dl: False
            mcmc (unweighted)
                    (416.8858712424856, 114) 72
            mcmc
                    (543.6668155389088, 115) 65
                    (-116.88680318312439, 114) 78
                    (64.10454898532073, 114) 70
            merge
                    (344.14747553149675, 50)
                    (395.5036614885022, 112)
                    (-144.23929897210388, 115)
            shrink
                    5

   directed: False overlap: False layered: True deg-corr: True dl: True
            mcmc (unweighted)
                    (-515.0632865530644, 114) 71
            mcmc
                    (-429.52721402339705, 112) 75
                    (-97.8538715008062, 110) 69
                    (61.234048124156786, 115) 73
            merge
                    (-750.3186792142322, 50)
                    (-289.89054259425507, 109)
                    (117.72715613340667, 115)
            shrink
                    5

   directed: False overlap: True layered: False deg-corr: False dl: False
            mcmc (unweighted)
                    (668.3427656597461, 1225) 785
                    (25.29309657374884, 1337) 773
            merge
                    (0.0, 50)
            shrink
                    5

   directed: False overlap: True layered: False deg-corr: False dl: True
            mcmc (unweighted)
                    (-454.0663525173459, 1223) 770
                    (-3.768593665481344, 1341) 772
            merge
                    (-302.1438895586243, 50)
            shrink
                    5

   directed: False overlap: True layered: False deg-corr: True dl: False
            mcmc (unweighted)
                    (715.484682404301, 1225) 765
                    (-34.61523962571205, 1341) 784
            merge
                    (34.77514206365362, 50)
            shrink
                    5

   directed: False overlap: True layered: False deg-corr: True dl: True
            mcmc (unweighted)
                    (51.108409710889575, 1226) 775
                    (-5.747478822563901, 1338) 771
            merge
                    (-288.7443154032152, 50)
            shrink
                    5

   directed: False overlap: True layered: covariates deg-corr: False dl: False
            mcmc (unweighted)
                    (30.51641860515133, 1198) 115
                    (-10.213872530535728, 1312) 115
            merge
                    (-279.91158404351484, 50)
            shrink
                    5

   directed: False overlap: True layered: covariates deg-corr: False dl: True
            mcmc (unweighted)
                    (-47.83809678456382, 1197) 115
                    (30.994946483584556, 1306) 115
            merge
                    (-1451.2844212406185, 50)
            shrink
                    5

   directed: False overlap: True layered: covariates deg-corr: True dl: False
            mcmc (unweighted)
                    (17.114695467713297, 1193) 115
                    (6.74992551722728, 1308) 115
            merge
                    (-141.66901107706147, 50)
            shrink
                    5

   directed: False overlap: True layered: covariates deg-corr: True dl: True
            mcmc (unweighted)
                    (18.28144731527671, 1186) 115
                    (-8.899892485955935, 1320) 115
            merge
                    (-2351.4411047024337, 50)
            shrink
                    5

   directed: False overlap: True layered: True deg-corr: False dl: False
            mcmc (unweighted)
                    (-30.414601069667068, 1192) 115
                    (23.38533942611877, 1317) 115
            merge
                    (-160.6578412585288, 50)
            shrink
                    5

   directed: False overlap: True layered: True deg-corr: False dl: True
            mcmc (unweighted)
                    (-15.9501258457887, 1185) 115
                    (43.8463674972257, 1312) 115
            merge
                    (-1100.5828201139577, 50)
            shrink
                    5

   directed: False overlap: True layered: True deg-corr: True dl: False
            mcmc (unweighted)
                    (38.98178711619759, 1197) 115
                    (18.742075941875243, 1319) 115
            merge
                    (-116.24155086445886, 50)
            shrink
                    5

   directed: False overlap: True layered: True deg-corr: True dl: True
            mcmc (unweighted)
                    (-0.7081085763732451, 1190) 115
                    (-31.565140167020655, 1307) 115
            merge
                    (-1350.8161986552966, 50)
            shrink
                    5

   directed: True overlap: False layered: False deg-corr: False
           Current bracket: (1, 60, 115) (2495.8290836444471, 3204.2321625907725, 3708.414497288833)
           Current bracket: (1, 26, 60) (2495.8290836444471, 2617.8774215566509, 3204.2321625907725)
           Current bracket: (1, 13, 26) (2495.8290836444471, 2363.1981701253062, 2617.8774215566509)
           Current bracket: (1, 13, 26) (2495.8290836444471, 2363.1981701253062, 2617.8774215566509)
           Bisect at B = 18 with S = 2453.367089312351
           Current bracket: (1, 13, 18) (2495.8290836444471, 2363.1981701253062, 2453.3670893123513)
           Bisect at B = 8 with S = 2300.994045668983
           Current bracket: (1, 8, 13) (2495.8290836444471, 2300.9940456689828, 2363.1981701253062)
           Bisect at B = 5 with S = 2285.052011050861
           Current bracket: (1, 5, 8) (2495.8290836444471, 2285.052011050861, 2300.9940456689828)
           Bisect at B = 3 with S = 2317.091178612431
           Current bracket: (3, 5, 8) (2317.0911786124307, 2285.052011050861, 2300.9940456689828)
           Bisect at B = 6 with S = 2284.866879995794
           Current bracket: (5, 6, 8) (2285.052011050861, 2284.8668799957941, 2300.9940456689828)
           Bisect at B = 7 with S = 2294.655787376431
           Current bracket: (5, 6, 7) (2285.052011050861, 2284.8668799957941, 2294.6557873764309)
           Bisect at B = 5 with S = 2285.052011050861
           Best result: B = 6, S = 2284.866879995794
   6 2284.86688
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 5) , dS: -210.777072594 2
           level 1 : replaced (5, 1) -> (5, 2) , dS: -2.66350792975 2
           level 2 : rejected replacement (2, 1) -> (2, 1) , dS: 0.0
           level 2 : rejected insert 2 , dS: 1.79175946923
           level 1 : skipping 2
           level 0 : replaced (115, 5) -> (115, 8) , dS: -5.91339284569 3
           level 1 : rejected replacement (8, 2) -> (8, 2) , dS: 0.0
           level 1 : rejected insert 8 , dS: 19.3741100228
           level 0 : skipping 8
   l: 0, N: 115, B: 8
   l: 1, N: 8, B: 2
   l: 2, N: 2, B: 1
   2276.47511028

   directed: True overlap: False layered: False deg-corr: True
           Current bracket: (1, 60, 115) (2318.8696624551153, 3051.6877329249223, 3337.5708321100378)
           Current bracket: (1, 26, 60) (2318.8696624551153, 2579.1905163091492, 3051.6877329249223)
           Current bracket: (1, 13, 26) (2318.8696624551153, 2292.6909131845441, 2579.1905163091492)
           Current bracket: (1, 13, 26) (2318.8696624551153, 2292.6909131845441, 2579.1905163091492)
           Bisect at B = 18 with S = 2424.807155254194
           Current bracket: (1, 13, 18) (2318.8696624551153, 2292.6909131845441, 2424.8071552541942)
           Bisect at B = 8 with S = 2196.047650916988
           Current bracket: (1, 8, 13) (2318.8696624551153, 2196.0476509169875, 2292.6909131845441)
           Bisect at B = 5 with S = 2165.076125813918
           Current bracket: (1, 5, 8) (2318.8696624551153, 2165.0761258139182, 2196.0476509169875)
           Bisect at B = 3 with S = 2178.869483800998
           Current bracket: (3, 5, 8) (2178.8694838009978, 2165.0761258139182, 2196.0476509169875)
           Bisect at B = 6 with S = 2170.665883964396
           Current bracket: (3, 5, 6) (2178.8694838009978, 2165.0761258139182, 2170.665883964396)
           Bisect at B = 4 with S = 2164.554037334525
           Best result: B = 4, S = 2164.554037334525
   4 2164.55403733
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 4) , dS: -158.455012563 2
           level 1 : replaced (4, 1) -> (4, 2) , dS: -0.230629487184 2
           level 2 : rejected replacement (2, 1) -> (2, 1) , dS: 0.0
           level 2 : rejected insert 2 , dS: 1.79175946923
           level 1 : skipping 2
           level 0 : rejected replacement (115, 4) -> (115, 5) , dS: 1.88971212096
   l: 0, N: 115, B: 4
   l: 1, N: 4, B: 2
   l: 2, N: 2, B: 1
   2160.1840204

   directed: True overlap: False layered: covariates deg-corr: False
           Current bracket: (1, 60, 115) (3876.323037275662, 4770.369390710639, 5076.3584970855154)
           Current bracket: (1, 26, 60) (3876.323037275662, 4390.4129518186537, 4770.369390710639)
           Current bracket: (1, 13, 26) (3876.323037275662, 4090.8678102257359, 4390.4129518186537)
           Current bracket: (1, 8, 13) (3876.323037275662, 3928.4123438440784, 4090.8678102257359)
           Current bracket: (1, 5, 8) (3876.323037275662, 3815.3301172970105, 3928.4123438440784)
           Current bracket: (1, 5, 8) (3876.323037275662, 3815.3301172970105, 3928.4123438440784)
           Bisect at B = 3 with S = 3770.207464051877
           Current bracket: (1, 3, 5) (3876.323037275662, 3770.2074640518772, 3815.3301172970105)
           Bisect at B = 2 with S = 3787.213321435399
           Current bracket: (2, 3, 5) (3787.2133214353989, 3770.2074640518772, 3815.3301172970105)
           Bisect at B = 4 with S = 3781.892332292212
           Current bracket: (2, 3, 4) (3787.2133214353989, 3770.2074640518772, 3781.8923322922124)
           Bisect at B = 2 with S = 3787.213321435399
           Best result: B = 3, S = 3770.207464051877
   3 3770.20746405
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 6) , dS: -201.16791632 2
           level 1 : rejected replacement (6, 1) -> (6, 1) , dS: 0.0
           level 1 : rejected insert 6 , dS: 12.7148161031
           level 0 : skipping 6
   l: 0, N: 115, B: 6
   l: 1, N: 6, B: 1
   3675.15512096

   directed: True overlap: False layered: covariates deg-corr: True
           Current bracket: (1, 60, 115) (3699.3636160863302, 4524.7650859192672, 4705.5148319067212)
           Current bracket: (1, 26, 60) (3699.3636160863302, 4280.5863383999085, 4524.7650859192672)
           Current bracket: (1, 13, 26) (3699.3636160863302, 4013.0275490563249, 4280.5863383999085)
           Current bracket: (1, 8, 13) (3699.3636160863302, 3825.7420416344448, 4013.0275490563249)
           Current bracket: (1, 5, 8) (3699.3636160863302, 3716.3154093581552, 3825.7420416344448)
           Current bracket: (1, 3, 5) (3699.3636160863302, 3636.3009795251073, 3716.3154093581552)
           Current bracket: (1, 3, 5) (3699.3636160863302, 3636.3009795251073, 3716.3154093581552)
           Bisect at B = 2 with S = 3638.452352760167
           Current bracket: (2, 3, 5) (3638.4523527601673, 3636.3009795251073, 3716.3154093581552)
           Bisect at B = 4 with S = 3673.832571602863
           Current bracket: (2, 3, 4) (3638.4523527601673, 3636.3009795251073, 3673.8325716028635)
           Bisect at B = 2 with S = 3638.452352760167
           Best result: B = 3, S = 3636.300979525107
   3 3636.30097953
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 5) , dS: -150.163186164 2
           level 1 : rejected replacement (5, 1) -> (5, 1) , dS: 0.0
           level 1 : rejected insert 5 , dS: 9.62377364973
           level 0 : skipping 5
   l: 0, N: 115, B: 5
   l: 1, N: 5, B: 1
   3549.20042992

   directed: True overlap: False layered: True deg-corr: False
           Current bracket: (1, 60, 115) (4149.3841336310379, 5284.3845874087538, 5806.3692226513531)
           Current bracket: (1, 26, 60) (4149.3841336310379, 4970.2156440906747, 5284.3845874087538)
           Current bracket: (1, 13, 26) (4149.3841336310379, 4665.913218701382, 4970.2156440906747)
           Current bracket: (1, 8, 13) (4149.3841336310379, 4448.1593590070361, 4665.913218701382)
           Current bracket: (1, 5, 8) (4149.3841336310379, 4259.415534978637, 4448.1593590070361)
           Current bracket: (1, 3, 5) (4149.3841336310379, 4156.7168613952836, 4259.415534978637)
           Current bracket: (1, 2, 3) (4149.3841336310379, 4108.7083204548444, 4156.7168613952836)
           Bisect at B = 1 with S = 4149.384133631038
           Best result: B = 2, S = 4108.708320454844
   2 4108.70832045
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : rejected replacement (115, 1) -> (115, 1) , dS: 0.0
   l: 0, N: 115, B: 1
   l: 1, N: 1, B: 1
   4151.68671872

   directed: True overlap: False layered: True deg-corr: True
           Current bracket: (1, 60, 115) (4079.1051067893313, 4596.3061341987523, 4808.8764024295706)
           Current bracket: (1, 26, 60) (4079.1051067893313, 4496.5820872632921, 4596.3061341987523)
           Current bracket: (1, 13, 26) (4079.1051067893313, 4370.6873281497483, 4496.5820872632921)
           Current bracket: (1, 8, 13) (4079.1051067893313, 4273.0526544098466, 4370.6873281497483)
           Current bracket: (1, 5, 8) (4079.1051067893313, 4185.3762872671505, 4273.0526544098466)
           Current bracket: (1, 3, 5) (4079.1051067893313, 4116.5895360878512, 4185.3762872671505)
           Current bracket: (1, 2, 3) (4079.1051067893313, 4113.7990467107265, 4116.5895360878512)
           Bisect at B = 1 with S = 4079.105106789331
           Best result: B = 1, S = 4079.105106789331
   1 4079.10510679
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : rejected replacement (115, 1) -> (115, 1) , dS: 0.0
   l: 0, N: 115, B: 1
   l: 1, N: 1, B: 1
   4081.40769188

   directed: True overlap: True layered: False deg-corr: False
           Current bracket: (1, 60, 115) (2495.8290836444471, 3243.3591987389864, 3708.4144972888325)
           Current bracket: (1, 26, 60) (2495.8290836444471, 2673.966224376918, 3243.3591987389864)
           Current bracket: (1, 13, 26) (2495.8290836444471, 2441.2363917608513, 2673.966224376918)
           Current bracket: (1, 13, 26) (2495.8290836444471, 2441.2363917608513, 2673.966224376918)
           Bisect at B = 18 with S = 2510.597188861897
           Current bracket: (1, 13, 18) (2495.8290836444471, 2441.2363917608513, 2510.5971888618969)
           Bisect at B = 8 with S = 2337.458764563424
           Current bracket: (1, 8, 13) (2495.8290836444471, 2337.4587645634238, 2441.2363917608513)
           Bisect at B = 5 with S = 2301.092707087893
           Current bracket: (1, 5, 8) (2495.8290836444471, 2301.0927070878934, 2337.4587645634238)
           Bisect at B = 3 with S = 2357.912801442382
           Current bracket: (3, 5, 8) (2357.9128014423818, 2301.0927070878934, 2337.4587645634238)
           Bisect at B = 6 with S = 2309.754290101913
           Current bracket: (3, 5, 6) (2357.9128014423818, 2301.0927070878934, 2309.7542901019133)
           Bisect at B = 4 with S = 2298.969652099869
           Best result: B = 4, S = 2298.969652099869
   4 2298.9696521
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 4) , dS: -167.727455603 2
           level 1 : rejected replacement (4, 1) -> (4, 1) , dS: 0.0
           level 1 : rejected insert 4 , dS: 6.73340189184
           level 0 : skipping 4
   l: 0, N: 115, B: 4
   l: 1, N: 4, B: 1
   2328.10162804

   directed: True overlap: True layered: False deg-corr: True
           Current bracket: (1, 60, 115) (2318.8696624551112, 3614.0932889394471, 3337.5708321100401)
           Current bracket: (1, 26, 60) (2318.8696624551112, 3110.1951905894075, 3614.0932889394471)
           Current bracket: (1, 13, 26) (2318.8696624551112, 2710.2613205743455, 3110.1951905894075)
           Current bracket: (1, 8, 13) (2318.8696624551112, 2528.3426363372578, 2710.2613205743455)
           Current bracket: (1, 5, 8) (2318.8696624551112, 2440.4089835968939, 2528.3426363372578)
           Current bracket: (1, 3, 5) (2318.8696624551112, 2353.166621031437, 2440.4089835968939)
           Current bracket: (1, 2, 3) (2318.8696624551112, 2264.6889302483723, 2353.166621031437)
           Bisect at B = 1 with S = 2318.869662455111
           Best result: B = 2, S = 2264.688930248372
   2 2264.68893025
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 2) , dS: -55.0106699939 2
           level 1 : rejected replacement (2, 1) -> (2, 1) , dS: 0.0
           level 1 : rejected insert 2 , dS: 1.79175946923
           level 0 : skipping 2
   l: 0, N: 115, B: 2
   l: 1, N: 2, B: 1
   2263.85899246

   directed: True overlap: True layered: covariates deg-corr: False
           Current bracket: (1, 60, 115) (3876.323037275662, 4795.4571396071951, 5076.3584970855154)
           Current bracket: (1, 26, 60) (3876.323037275662, 4480.189057603242, 4795.4571396071951)
           Current bracket: (1, 13, 26) (3876.323037275662, 4232.1382878497698, 4480.189057603242)
           Current bracket: (1, 8, 13) (3876.323037275662, 4069.039478169384, 4232.1382878497698)
           Current bracket: (1, 5, 8) (3876.323037275662, 3951.8135764459576, 4069.039478169384)
           Current bracket: (1, 3, 5) (3876.323037275662, 3878.0269892829961, 3951.8135764459576)
           Current bracket: (1, 2, 3) (3876.323037275662, 3867.7472612920978, 3878.0269892829961)
           Bisect at B = 1 with S = 3876.323037275662
           Best result: B = 2, S = 3867.747261292098
   2 3867.74726129
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 5) , dS: -178.717049511 2
           level 1 : rejected replacement (5, 1) -> (5, 1) , dS: 0.0
           level 1 : rejected insert 5 , dS: 9.62377364973
           level 0 : skipping 5
   l: 0, N: 115, B: 5
   l: 1, N: 5, B: 1
   3697.60598777

   directed: True overlap: True layered: covariates deg-corr: True
           Current bracket: (1, 60, 115) (3699.3636160863266, 5167.062731660124, 4705.514831906723)
           Current bracket: (1, 26, 60) (3699.3636160863266, 4845.8008598492379, 5167.062731660124)
           Current bracket: (1, 13, 26) (3699.3636160863266, 4614.9766720583375, 4845.8008598492379)
           Current bracket: (1, 8, 13) (3699.3636160863266, 4368.643894604631, 4614.9766720583375)
           Current bracket: (1, 5, 8) (3699.3636160863266, 4154.7195000910806, 4368.643894604631)
           Current bracket: (1, 3, 5) (3699.3636160863266, 3929.3119299306791, 4154.7195000910806)
           Current bracket: (1, 2, 3) (3699.3636160863266, 3780.6025531339842, 3929.3119299306791)
           Bisect at B = 1 with S = 3699.363616086327
           Best result: B = 1, S = 3699.363616086327
   1 3699.36361609
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 5) , dS: -135.304469399 2
           level 1 : rejected replacement (5, 1) -> (5, 1) , dS: 0.0
           level 1 : rejected insert 5 , dS: 9.62377364973
           level 0 : skipping 5
   l: 0, N: 115, B: 5
   l: 1, N: 5, B: 1
   3564.05914669

   directed: True overlap: True layered: True deg-corr: False
           Current bracket: (1, 60, 115) (4149.3841336310379, 5634.8056324677736, 5806.3692226513531)
           Current bracket: (1, 26, 60) (4149.3841336310379, 5253.8064287341322, 5634.8056324677736)
           Current bracket: (1, 13, 26) (4149.3841336310379, 4830.7710230477469, 5253.8064287341322)
           Current bracket: (1, 8, 13) (4149.3841336310379, 4295.682486401598, 4830.7710230477469)
           Current bracket: (1, 5, 8) (4149.3841336310379, 4052.3861265593323, 4295.682486401598)
           Current bracket: (1, 5, 8) (4149.3841336310379, 4052.3861265593323, 4295.682486401598)
           Bisect at B = 3 with S = 3833.118315845553
           Current bracket: (1, 3, 5) (4149.3841336310379, 3833.1183158455533, 4052.3861265593323)
           Bisect at B = 2 with S = 4039.690876848615
           Current bracket: (2, 3, 5) (4039.6908768486151, 3833.1183158455533, 4052.3861265593323)
           Bisect at B = 4 with S = 4036.591921918258
           Current bracket: (2, 3, 4) (4039.6908768486151, 3833.1183158455533, 4036.5919219182579)
           Bisect at B = 2 with S = 4039.690876848615
           Best result: B = 3, S = 3833.118315845553
   3 3833.11831585
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : rejected replacement (115, 1) -> (115, 1) , dS: 0.0
   l: 0, N: 115, B: 1
   l: 1, N: 1, B: 1
   4151.68671872

   directed: True overlap: True layered: True deg-corr: True
           Current bracket: (1, 60, 115) (4079.1051067893313, 5197.2037808565174, 4808.8764024295706)
           Current bracket: (1, 26, 60) (4079.1051067893313, 4947.0264093741907, 5197.2037808565174)
           Current bracket: (1, 13, 26) (4079.1051067893313, 4812.7634310564736, 4947.0264093741907)
           Current bracket: (1, 8, 13) (4079.1051067893313, 4602.6320859157631, 4812.7634310564736)
           Current bracket: (1, 5, 8) (4079.1051067893313, 4252.8246786173822, 4602.6320859157631)
           Current bracket: (1, 3, 5) (4079.1051067893313, 4049.5683588147385, 4252.8246786173822)
           Current bracket: (1, 3, 5) (4079.1051067893313, 4049.5683588147385, 4252.8246786173822)
           Bisect at B = 2 with S = 4074.521377773606
           Current bracket: (2, 3, 5) (4074.5213777736058, 4049.5683588147385, 4252.8246786173822)
           Bisect at B = 4 with S = 4118.559112781375
           Current bracket: (2, 3, 4) (4074.5213777736058, 4049.5683588147385, 4118.559112781375)
           Bisect at B = 2 with S = 4074.521377773606
           Best result: B = 3, S = 4049.568358814739
   3 4049.56835881
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : rejected replacement (115, 1) -> (115, 1) , dS: 0.0
   l: 0, N: 115, B: 1
   l: 1, N: 1, B: 1
   4081.40769188

   directed: False overlap: False layered: False deg-corr: False
           Current bracket: (1, 60, 115) (2070.9298619612, 2742.0286346974412, 3302.1523929668851)
           Current bracket: (1, 26, 60) (2070.9298619612, 2163.2011864326296, 2742.0286346974412)
           Current bracket: (1, 13, 26) (2070.9298619612, 1884.6261937918878, 2163.2011864326296)
           Current bracket: (1, 13, 26) (2070.9298619612, 1884.6261937918878, 2163.2011864326296)
           Bisect at B = 18 with S = 1991.644976112192
           Current bracket: (1, 13, 18) (2070.9298619612, 1884.6261937918878, 1991.6449761121924)
           Bisect at B = 8 with S = 1824.47791693959
           Current bracket: (1, 8, 13) (2070.9298619612, 1824.4779169395899, 1884.6261937918878)
           Bisect at B = 5 with S = 1832.933393330268
           Current bracket: (5, 8, 13) (1832.933393330268, 1824.4779169395899, 1884.6261937918878)
           Bisect at B = 10 with S = 1833.549647859825
           Current bracket: (5, 8, 10) (1832.933393330268, 1824.4779169395899, 1833.5496478598247)
           Bisect at B = 6 with S = 1824.864187113754
           Current bracket: (6, 8, 10) (1824.8641871137543, 1824.4779169395899, 1833.5496478598247)
           Bisect at B = 7 with S = 1821.898539100332
           Current bracket: (6, 7, 8) (1824.8641871137543, 1821.8985391003318, 1824.4779169395899)
           Bisect at B = 6 with S = 1824.864187113754
           Best result: B = 7, S = 1821.898539100332
   7 1821.8985391
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 10) , dS: -247.539097789 2
           level 1 : replaced (10, 1) -> (10, 2) , dS: -8.57441720825 2
           level 2 : rejected replacement (2, 1) -> (2, 1) , dS: 0.0
           level 2 : rejected insert 2 , dS: 1.79175946923
           level 1 : skipping 2
           level 0 : rejected replacement (115, 10) -> (115, 9) , dS: 1.04063514461
   l: 0, N: 115, B: 10
   l: 1, N: 10, B: 2
   l: 2, N: 2, B: 1
   1814.81634696

   directed: False overlap: False layered: False deg-corr: True
           Current bracket: (1, 60, 115) (2039.438415802616, 2687.4979970801901, 3059.7072673586936)
           Current bracket: (1, 26, 60) (2039.438415802616, 2210.1885588389487, 2687.4979970801901)
           Current bracket: (1, 13, 26) (2039.438415802616, 1956.4109872858012, 2210.1885588389487)
           Current bracket: (1, 13, 26) (2039.438415802616, 1956.4109872858012, 2210.1885588389487)
           Bisect at B = 18 with S = 2056.248563407142
           Current bracket: (1, 13, 18) (2039.438415802616, 1956.4109872858012, 2056.2485634071418)
           Bisect at B = 8 with S = 1887.564759326295
           Current bracket: (1, 8, 13) (2039.438415802616, 1887.5647593262947, 1956.4109872858012)
           Bisect at B = 5 with S = 1878.123470790204
           Current bracket: (1, 5, 8) (2039.438415802616, 1878.1234707902045, 1887.5647593262947)
           Bisect at B = 3 with S = 1900.997091163706
           Current bracket: (3, 5, 8) (1900.9970911637056, 1878.1234707902045, 1887.5647593262947)
           Bisect at B = 6 with S = 1873.870847992987
           Current bracket: (5, 6, 8) (1878.1234707902045, 1873.8708479929869, 1887.5647593262947)
           Bisect at B = 7 with S = 1878.609827955009
           Current bracket: (5, 6, 7) (1878.1234707902045, 1873.8708479929869, 1878.6098279550095)
           Bisect at B = 5 with S = 1878.123470790204
           Best result: B = 6, S = 1873.870847992987
   6 1873.87084799
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 5) , dS: -164.055568622 2
           level 1 : rejected replacement (5, 1) -> (5, 1) , dS: 0.0
           level 1 : rejected insert 5 , dS: 9.62377364973
           level 0 : skipping 5
   l: 0, N: 115, B: 5
   l: 1, N: 5, B: 1
   1875.38284718

   directed: False overlap: False layered: covariates deg-corr: False
           Current bracket: (1, 60, 115) (3451.4238155924145, 4295.9957091835586, 4658.1574421850237)
           Current bracket: (1, 26, 60) (3451.4238155924145, 3824.5164704145741, 4295.9957091835586)
           Current bracket: (1, 13, 26) (3451.4238155924145, 3481.1896063728163, 3824.5164704145741)
           Current bracket: (1, 8, 13) (3451.4238155924145, 3350.3690166608876, 3481.1896063728163)
           Current bracket: (1, 8, 13) (3451.4238155924145, 3350.3690166608876, 3481.1896063728163)
           Bisect at B = 5 with S = 3305.96690163836
           Current bracket: (1, 5, 8) (3451.4238155924145, 3305.9669016383605, 3350.3690166608876)
           Bisect at B = 3 with S = 3305.958349260824
           Current bracket: (1, 3, 5) (3451.4238155924145, 3305.958349260824, 3305.9669016383605)
           Bisect at B = 2 with S = 3344.118688610963
           Current bracket: (2, 3, 5) (3344.1186886109631, 3305.958349260824, 3305.9669016383605)
           Bisect at B = 4 with S = 3290.938676531715
           Current bracket: (3, 4, 5) (3305.958349260824, 3290.9386765317154, 3305.9669016383605)
           Bisect at B = 3 with S = 3305.958349260824
           Best result: B = 4, S = 3290.938676531715
   4 3290.93867653
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 9) , dS: -236.856154187 2
           level 1 : rejected replacement (9, 1) -> (9, 1) , dS: 0.0
           level 1 : rejected insert 9 , dS: 22.9004705474
           level 0 : skipping 9
   l: 0, N: 115, B: 9
   l: 1, N: 9, B: 1
   3214.56766141

   directed: False overlap: False layered: covariates deg-corr: True
           Current bracket: (1, 60, 115) (3419.9323694338309, 4223.9918797736118, 4415.7123165768317)
           Current bracket: (1, 26, 60) (3419.9323694338309, 3875.1607535720391, 4223.9918797736118)
           Current bracket: (1, 13, 26) (3419.9323694338309, 3553.4454991243847, 3875.1607535720391)
           Current bracket: (1, 8, 13) (3419.9323694338309, 3406.1491969341596, 3553.4454991243847)
           Current bracket: (1, 8, 13) (3419.9323694338309, 3406.1491969341596, 3553.4454991243847)
           Bisect at B = 5 with S = 3330.744584199459
           Current bracket: (1, 5, 8) (3419.9323694338309, 3330.7445841994586, 3406.1491969341596)
           Bisect at B = 3 with S = 3317.96165911753
           Current bracket: (1, 3, 5) (3419.9323694338309, 3317.9616591175295, 3330.7445841994586)
           Bisect at B = 2 with S = 3338.6229117634
           Current bracket: (2, 3, 5) (3338.6229117633998, 3317.9616591175295, 3330.7445841994586)
           Bisect at B = 4 with S = 3318.361316320577
           Current bracket: (2, 3, 4) (3338.6229117633998, 3317.9616591175295, 3318.3613163205769)
           Bisect at B = 2 with S = 3338.6229117634
           Best result: B = 3, S = 3317.96165911753
   3 3317.96165912
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 6) , dS: -152.814489428 2
           level 1 : rejected replacement (6, 1) -> (6, 1) , dS: 0.0
           level 1 : rejected insert 6 , dS: 12.7148161031
           level 0 : skipping 6
   l: 0, N: 115, B: 6
   l: 1, N: 6, B: 1
   3267.11788001

   directed: False overlap: False layered: True deg-corr: False
           Current bracket: (1, 60, 115) (3724.4849119477913, 4846.7272154433076, 5388.1681677508614)
           Current bracket: (1, 26, 60) (3724.4849119477913, 4415.9298484742922, 4846.7272154433076)
           Current bracket: (1, 13, 26) (3724.4849119477913, 4043.547617400895, 4415.9298484742922)
           Current bracket: (1, 8, 13) (3724.4849119477913, 3831.4149302559326, 4043.547617400895)
           Current bracket: (1, 5, 8) (3724.4849119477913, 3716.0144978844455, 3831.4149302559326)
           Current bracket: (1, 5, 8) (3724.4849119477913, 3716.0144978844455, 3831.4149302559326)
           Bisect at B = 3 with S = 3664.900838079229
           Current bracket: (1, 3, 5) (3724.4849119477913, 3664.9008380792288, 3716.0144978844455)
           Bisect at B = 2 with S = 3657.203439841218
           Current bracket: (1, 2, 3) (3724.4849119477913, 3657.2034398412184, 3664.9008380792288)
           Bisect at B = 1 with S = 3724.484911947791
           Best result: B = 2, S = 3657.203439841218
   2 3657.20343984
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : replaced (115, 1) -> (115, 2) , dS: -52.0783689941 2
           level 1 : rejected replacement (2, 1) -> (2, 1) , dS: 0.0
           level 1 : rejected insert 2 , dS: 5.79909265446
           level 0 : skipping 2
   l: 0, N: 115, B: 2
   l: 1, N: 2, B: 1
   3674.70912805

   directed: False overlap: False layered: True deg-corr: True
           Current bracket: (1, 60, 115) (3690.8748461232653, 4390.7235884784623, 4494.4233125699411)
           Current bracket: (1, 26, 60) (3690.8748461232653, 4212.4807663850643, 4390.7235884784623)
           Current bracket: (1, 13, 26) (3690.8748461232653, 4029.1609658456218, 4212.4807663850643)
           Current bracket: (1, 8, 13) (3690.8748461232653, 3913.6020693472092, 4029.1609658456218)
           Current bracket: (1, 5, 8) (3690.8748461232653, 3834.8957858009617, 3913.6020693472092)
           Current bracket: (1, 3, 5) (3690.8748461232653, 3727.2813995894244, 3834.8957858009617)
           Current bracket: (1, 2, 3) (3690.8748461232653, 3690.288221840598, 3727.2813995894244)
           Bisect at B = 1 with S = 3690.874846123265
           Best result: B = 2, S = 3690.288221840598
   2 3690.28822184
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : rejected replacement (115, 1) -> (115, 1) , dS: 0.0
   l: 0, N: 115, B: 1
   l: 1, N: 1, B: 1
   3693.17743122

   directed: False overlap: True layered: False deg-corr: False
           Current bracket: (1, 60, 115) (2070.9298619612, 2773.8103951391176, 3302.1523929668847)
           Current bracket: (1, 26, 60) (2070.9298619612, 2197.5369172261844, 2773.8103951391176)
           Current bracket: (1, 13, 26) (2070.9298619612, 1888.693285109395, 2197.5369172261844)
           Current bracket: (1, 13, 26) (2070.9298619612, 1888.693285109395, 2197.5369172261844)
           Bisect at B = 18 with S = 2008.358125927789
           Current bracket: (1, 13, 18) (2070.9298619612, 1888.693285109395, 2008.3581259277889)
           Bisect at B = 8 with S = 1821.812463004034
           Current bracket: (1, 8, 13) (2070.9298619612, 1821.8124630040336, 1888.693285109395)
           Bisect at B = 5 with S = 1853.367979967209
           Current bracket: (5, 8, 13) (1853.3679799672091, 1821.8124630040336, 1888.693285109395)
           Bisect at B = 10 with S = 1823.97846127917
           Current bracket: (5, 8, 10) (1853.3679799672091, 1821.8124630040336, 1823.9784612791696)
           Bisect at B = 6 with S = 1828.525381763233
           Current bracket: (6, 8, 10) (1828.5253817632329, 1821.8124630040336, 1823.9784612791696)
           Bisect at B = 7 with S = 1823.245540811884
           Current bracket: (7, 8, 10) (1823.2455408118842, 1821.8124630040336, 1823.9784612791696)
           Bisect at B = 9 with S = 1822.086944968345
           Current bracket: (7, 8, 9) (1823.2455408118842, 1821.8124630040336, 1822.0869449683451)
           Bisect at B = 7 with S = 1823.245540811884
           Best result: B = 8, S = 1821.812463004034
   8 1821.812463
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 8) , dS: -232.661496277 2
           level 1 : replaced (8, 1) -> (8, 2) , dS: -2.95970180508 2
           level 2 : rejected replacement (2, 1) -> (2, 1) , dS: 0.0
           level 2 : rejected insert 2 , dS: 1.79175946923
           level 1 : skipping 2
           level 0 : replaced (115, 8) -> (115, 8) , dS: -4.93282634306 3
           level 1 : rejected replacement (8, 2) -> (8, 2) , dS: 0.0
           level 1 : rejected insert 8 , dS: 19.3741100228
           level 0 : skipping 8
   l: 0, N: 115, B: 8
   l: 1, N: 8, B: 2
   l: 2, N: 2, B: 1
   1830.37583754

   directed: False overlap: True layered: False deg-corr: True
           Current bracket: (1, 60, 115) (2039.438415802616, 2791.218655464937, 3059.7072673586931)
           Current bracket: (1, 26, 60) (2039.438415802616, 2596.1318468202285, 2791.218655464937)
           Current bracket: (1, 13, 26) (2039.438415802616, 2249.7395717949853, 2596.1318468202285)
           Current bracket: (1, 8, 13) (2039.438415802616, 2192.5567166861247, 2249.7395717949853)
           Current bracket: (1, 5, 8) (2039.438415802616, 2155.8270672199023, 2192.5567166861247)
           Current bracket: (1, 3, 5) (2039.438415802616, 2097.064562904091, 2155.8270672199023)
           Current bracket: (1, 2, 3) (2039.438415802616, 2067.0267866042686, 2097.064562904091)
           Bisect at B = 1 with S = 2039.438415802616
           Best result: B = 1, S = 2039.438415802616
   1 2039.4384158
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : rejected replacement (115, 1) -> (115, 1) , dS: 0.0
   l: 0, N: 115, B: 1
   l: 1, N: 1, B: 1
   2039.4384158

   directed: False overlap: True layered: covariates deg-corr: False
           Current bracket: (1, 60, 115) (3451.4238155924145, 4331.561281582568, 4658.1574421850237)
           Current bracket: (1, 26, 60) (3451.4238155924145, 3901.5868258047312, 4331.561281582568)
           Current bracket: (1, 13, 26) (3451.4238155924145, 3673.3309590473609, 3901.5868258047312)
           Current bracket: (1, 8, 13) (3451.4238155924145, 3552.8365768937738, 3673.3309590473609)
           Current bracket: (1, 5, 8) (3451.4238155924145, 3445.7715685052244, 3552.8365768937738)
           Current bracket: (1, 5, 8) (3451.4238155924145, 3445.7715685052244, 3552.8365768937738)
           Bisect at B = 3 with S = 3378.659618732226
           Current bracket: (1, 3, 5) (3451.4238155924145, 3378.6596187322261, 3445.7715685052244)
           Bisect at B = 2 with S = 3384.190161466892
           Current bracket: (2, 3, 5) (3384.1901614668918, 3378.6596187322261, 3445.7715685052244)
           Bisect at B = 4 with S = 3395.991818320718
           Current bracket: (2, 3, 4) (3384.1901614668918, 3378.6596187322261, 3395.9918183207178)
           Bisect at B = 2 with S = 3384.190161466892
           Best result: B = 3, S = 3378.659618732226
   3 3378.65961873
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 8) , dS: -236.467003238 2
           level 1 : rejected replacement (8, 1) -> (8, 1) , dS: 0.0
           level 1 : rejected insert 8 , dS: 19.3741100228
           level 0 : skipping 8
   l: 0, N: 115, B: 8
   l: 1, N: 8, B: 1
   3214.95681235

   directed: False overlap: True layered: covariates deg-corr: True
           Current bracket: (1, 60, 115) (3419.9323694338309, 4273.5270137403832, 4415.7123165768317)
           Current bracket: (1, 26, 60) (3419.9323694338309, 4003.9207484397484, 4273.5270137403832)
           Current bracket: (1, 13, 26) (3419.9323694338309, 3747.5119493111124, 4003.9207484397484)
           Current bracket: (1, 8, 13) (3419.9323694338309, 3587.8678820583082, 3747.5119493111124)
           Current bracket: (1, 5, 8) (3419.9323694338309, 3527.8285413311705, 3587.8678820583082)
           Current bracket: (1, 3, 5) (3419.9323694338309, 3450.2450966889314, 3527.8285413311705)
           Current bracket: (1, 2, 3) (3419.9323694338309, 3388.2681843019941, 3450.2450966889314)
           Bisect at B = 1 with S = 3419.932369433831
           Best result: B = 2, S = 3388.268184301994
   2 3388.2681843
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 0.0
           level 0 : replaced (115, 1) -> (115, 6) , dS: -164.980561445 2
           level 1 : rejected replacement (6, 1) -> (6, 1) , dS: 0.0
           level 1 : rejected insert 6 , dS: 12.7148161031
           level 0 : skipping 6
   l: 0, N: 115, B: 6
   l: 1, N: 6, B: 1
   3254.95180799

   directed: False overlap: True layered: True deg-corr: False
           Current bracket: (1, 60, 115) (3724.4849119477913, 5222.4917279727715, 5388.1681677508614)
           Current bracket: (1, 26, 60) (3724.4849119477913, 4745.7356091147667, 5222.4917279727715)
           Current bracket: (1, 13, 26) (3724.4849119477913, 4150.3055871227125, 4745.7356091147667)
           Current bracket: (1, 8, 13) (3724.4849119477913, 3778.92601352782, 4150.3055871227125)
           Current bracket: (1, 5, 8) (3724.4849119477913, 3587.5603537522038, 3778.92601352782)
           Current bracket: (1, 5, 8) (3724.4849119477913, 3587.5603537522038, 3778.92601352782)
           Bisect at B = 3 with S = 3529.792099317828
           Current bracket: (1, 3, 5) (3724.4849119477913, 3529.7920993178282, 3587.5603537522038)
           Bisect at B = 2 with S = 3545.492577998865
           Current bracket: (2, 3, 5) (3545.4925779988648, 3529.7920993178282, 3587.5603537522038)
           Bisect at B = 4 with S = 3552.535413509081
           Current bracket: (2, 3, 4) (3545.4925779988648, 3529.7920993178282, 3552.5354135090811)
           Bisect at B = 2 with S = 3545.492577998865
           Best result: B = 3, S = 3529.792099317828
   3 3529.79209932
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : replaced (115, 1) -> (115, 2) , dS: -52.0783689941 2
           level 1 : rejected replacement (2, 1) -> (2, 1) , dS: 0.0
           level 1 : rejected insert 2 , dS: 5.79909265446
           level 0 : skipping 2
   l: 0, N: 115, B: 2
   l: 1, N: 2, B: 1
   3674.70912805

   directed: False overlap: True layered: True deg-corr: True
           Current bracket: (1, 60, 115) (3690.8748461232653, 4757.8922348797005, 4494.4233125699411)
           Current bracket: (1, 26, 60) (3690.8748461232653, 4578.68710077166, 4757.8922348797005)
           Current bracket: (1, 13, 26) (3690.8748461232653, 4249.4783928541656, 4578.68710077166)
           Current bracket: (1, 8, 13) (3690.8748461232653, 3986.881326217177, 4249.4783928541656)
           Current bracket: (1, 5, 8) (3690.8748461232653, 3856.2684964836299, 3986.881326217177)
           Current bracket: (1, 3, 5) (3690.8748461232653, 3699.7810945700962, 3856.2684964836299)
           Current bracket: (1, 2, 3) (3690.8748461232653, 3641.6685995295234, 3699.7810945700962)
           Bisect at B = 1 with S = 3690.874846123265
           Best result: B = 2, S = 3641.668599529523
   2 3641.66859953
           level 1 : rejected replacement (1, 1) -> (1, 1) , dS: 0.0
           level 1 : rejected insert 1 , dS: 2.30258509299
           level 0 : rejected replacement (115, 1) -> (115, 1) , dS: 0.0
   l: 0, N: 115, B: 1
   l: 1, N: 1, B: 1
   3693.17743122
