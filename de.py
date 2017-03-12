


import sys




bdir1 = sys.argv[1]

bdir2= sys.argv[2]


foldchange = float(sys.argv[3])
bsupport = float(sys.argv[4])


mode = sys.argv[5]


prefix = "perm_"
nboot = 100


metabolic = ["ko01200","ko01210","ko01212","ko01230","ko01220","ko00010","ko00020","ko00030","ko00040","ko00051","ko00052","ko00053","ko00500","ko00520","ko00620","ko00630","ko00640","ko00650","ko00660","ko00562","ko00190","ko00195","ko00196","ko00710","ko00720","ko00680","ko00910","ko00920","ko00061","ko00062","ko00071","ko00072","ko00073","ko00100","ko00120","ko00121","ko00140","ko00561","ko00564","ko00565","ko00600","ko00590","ko00591","ko00592","ko01040","ko00230","ko00240","ko00250","ko00260","ko00270","ko00280","ko00290","ko00300","ko00310","ko00220","ko00330","ko00340","ko00350","ko00360","ko00380","ko00400","ko00410","ko00430","ko00440","ko00450","ko00460","ko00471","ko00472","ko00473","ko00480","ko00510","ko00513","ko00512","ko00515","ko00514","ko00532","ko00534","ko00533","ko00531","ko00563","ko00601","ko00603","ko00604","ko00540","ko00550","ko00511","ko00730","ko00740","ko00750","ko00760","ko00770","ko00780","ko00785","ko00790","ko00670","ko00830","ko00860","ko00130","ko00900","ko00902","ko00909","ko00904","ko00906","ko00905","ko00981","ko00908","ko00903","ko00281","ko01052","ko00522","ko01051","ko01059","ko01056","ko01057","ko00253","ko00523","ko01054","ko01053","ko01055","ko00940","ko00945","ko00941","ko00944","ko00942","ko00943","ko00901","ko00403","ko00950","ko00960","ko01058","ko00232","ko00965","ko00966","ko00402","ko00311","ko00332","ko00261","ko00331","ko00521","ko00524","ko00525","ko00231","ko00401","ko00404","ko00254","ko00362","ko00627","ko00364","ko00625","ko00361","ko00623","ko00622","ko00633","ko00642","ko00643","ko00791","ko00930","ko00351","ko00363","ko00621","ko00626","ko00624","ko00365","ko00984","ko00980","ko00982","ko00983","ko00194","ko01004","ko01007","ko01003","ko00535","ko00536","ko00537","ko01005","ko01006","ko01008","ko01000","ko01001","ko01009","ko01002","ko00199"]


def dictdict(bdir):
    ddict = {}
    for xx in range(1, nboot + 1):
        with open(bdir + "/perm_" + str(xx)) as f:
            a = f.readlines()
        a = map(lambda x: x.strip().split(), a)
        adict = {}
        for x, y in a:
            adict[x] = float(y)
        bdict = {}
        for x in metabolic:
            bdict[x] = adict.get(x, 0)
        ddict["perm_" + str(xx)] = bdict
    return ddict

exp1 = dictdict(bdir1)
exp2 = dictdict(bdir2)


metadict1 = {}
for x in metabolic:
    vals = []
    for i in range(1, nboot + 1):
        vals.append(exp1["perm_" + str(i)][x])
    metadict1[x] = vals

metadict2 = {}
for x in metabolic:
    vals = []
    for i in range(1, nboot + 1):
        vals.append(exp2["perm_" + str(i)][x])
    metadict2[x] = vals





for xx in metabolic:
    sample1 = metadict1[xx]
    sample2 = metadict2[xx]
    oll = 0
    fchange_good = 0
    if mode == "match":
        for x, y in zip(sample1, sample2):
            if x > foldchange * y:
                fchange_good += 1
            oll += 1
    elif mode == "all":
        for x in sample1:
            for y in sample2:
                if x > foldchange * y:
                    fchange_good += 1
                oll += 1
    if fchange_good * 1.0 / oll >= bsupport:
        print xx, max(sample1), max(sample2), fchange_good * 1.0 / oll





