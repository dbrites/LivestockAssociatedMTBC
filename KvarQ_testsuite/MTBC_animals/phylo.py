VERSION = '0.8'
from kvarq.genes import COMPATIBILITY as GENES_COMPATIBILITY

from kvarq.log import lo
from kvarq.genes import Reference, Test, SNP, Genotype, Testsuite

from _util import ancestor

class PhyloTestsuite(Testsuite):

    def __str__(self):
        return 'MTBC animal lineage SNPs'

    def score_SNPs(self, genotypes, coverages):
        ''' returns a dictionary ``{genotype:scores, ...}`` where ``scores`` is an
            array where every SNP for the given genotype is represented with a
            boolean value indicating whether the SNP is covered "sufficiently" 
            to be considered positive '''

        ret = {}

        for test in self.tests:

            coverage = coverages[test]
            genotype = test.genotype

            if genotype in genotypes:
                if test.template.validate(coverage):
                    ret.setdefault(genotype, []).append(True)
                else:
                    ret.setdefault(genotype, []).append(False)

        return ret

    def _analyse(self, coverages):
        mls = []

        #TODO choose criteria dynamically

        for ml, xs in self.score_SNPs(Lineage.roots, coverages).items():
            lo.debug(str(ml)+' : '+str(xs))

            if sum(xs)>1:
                # we need at least two positive SNPs
                mls.append(ml.name)

                if 0 in xs:
                    # flag if one of the SNPs is not found
                    # mls[-1] += ' (?)'
                    pass

                if ml.children:
                    sls = []

                    # co-mutants complicate our life somewhat
                    slsc = self.score_SNPs(ml.children, coverages)
                    slsc_byname = {}
                    slsc_comutants = {}
                    for sl, xs_ in slsc.items():
                        slsc_byname.setdefault(sl.name, []).extend(xs_)
                        if sl.comutant:
                            slsc_comutants.setdefault(sl.name, []).extend([sl.comutant] * sum(xs_))

                    for slname, xs_ in slsc_byname.items():
                        comutants = ''.join(slsc_comutants.get(slname, []))
                        lo.debug('sublineage '+slname+' : '+str(xs_)+' comutants '+comutants)
                        if sum(xs_)>1:
                            sls.append(slname)
                            if comutants:
                                sls[-1] += '_' + comutants

#                            if 0 in xs_: # does not make sense when using comutants
#                                sls[-1] += ' (?)'

                    if sls:
                        mls[-1] += '/' + '-'.join(sls)

        depths = sorted([coverage.mean(include_margins=False)
                for coverage in coverages.values()])
        remark = ''

        if depths[len(depths)/2] < 10:
            remark += ' -- low coverage (median below 10x)'

        mixed = sum([coverage.mixed() for coverage in coverages.values()])
        if mixed:
            remark += ' -- mixed coverage'

        if not mls:
            return '?' + remark

        return ' // '.join(mls) + remark


class Lineage(Genotype):

    roots = []

    def __init__(self, name, parent=None, comutant=None):
        super(Lineage, self).__init__(name)
        self.name = name
        self.parent = parent
        self.comutant = comutant
        self.children = []

        if parent:
            parent.children.append(self)
            
        else:
            Lineage.roots.append(self)


zwyer21 = Reference('Our reference')
comas09  = Reference('PLoS ONE 2009 - Comas (monomorphic)')
stucki12 = Reference('Stucki et al. PLoS ONE 2012')
ngabonziza20 = Reference('Ngabonziza et al. Nature communications 2020')

orygis         = Lineage('La3')
bovis         = Lineage('La1')
caprae  = Lineage('La2')
unk4_unk5_eu2         = Lineage('La1.7', bovis)
eu1_unk6_unk7        = Lineage('La1.8', bovis)
unk2         = Lineage('La1.2', bovis)
af2         = Lineage('La1.3', bovis)
af1         = Lineage('La1.6', bovis)
PZA_sus   = Lineage('La1.1', bovis)
unk3   = Lineage('La1.4', bovis)
eu1   = Lineage('La1.8.1', bovis)
unk4   = Lineage('unk4', bovis)
eu2   = Lineage('La1.7.1', bovis)
unk7   = Lineage('La1.8.2', bovis)
unk5   = Lineage('unk5', bovis)
unk6   = Lineage('unk6', bovis)
unk9   = Lineage('La1.5', bovis)
unk2_BCG = Lineage('La1.2_BCG', bovis)
lineage1         = Lineage('L1')
lineage2         = Lineage('L2')
lineage3         = Lineage('L3')
lineage4         = Lineage('L4')
lineage5         = Lineage('L5')
lineage6         = Lineage('L6')
lineage7         = Lineage('L7')
lineage8         = Lineage('L8')
lineage9         = Lineage('L9')



mainlineages_SNPs = [Test(SNP(genome=ancestor,pos=2475888,base='G'),bovis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=3877256,base='T'),bovis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=2612300,base='T'),bovis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=2912516,base='A'),bovis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=3830348,base='C'),orygis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=3396569,base='A'),orygis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=3577779,base='A'),orygis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=3770449,base='A'),orygis,zwyer21),
                     Test(SNP(genome=ancestor,pos=327232,base='C'),orygis,zwyer21), 
                     Test(SNP(genome=ancestor,pos=2382830,base='A'),caprae,zwyer21), 
                     Test(SNP(genome=ancestor,pos=2194204,base='G'),caprae,zwyer21), 
                     Test(SNP(genome=ancestor,pos=533328,base='A'),caprae,zwyer21), 
                     Test(SNP(genome=ancestor,pos=2701538,base='A'),caprae,zwyer21),
                     Test(SNP(genome=ancestor,pos=2731555,base='C'),caprae,zwyer21),
Test(SNP(genome=ancestor,pos=3920109,base='T'),lineage1,stucki12),
    Test(SNP(genome=ancestor,pos=3597682,base='T'),lineage1,stucki12),
    Test(SNP(genome=ancestor,pos=1590555,base='T'),lineage1,stucki12),
    Test(SNP(genome=ancestor,pos=1834177,base='C'),lineage2,stucki12),
    Test(SNP(genome=ancestor,pos=3304966,base='A'),lineage2,stucki12),
    Test(SNP(genome=ancestor,pos=2711722,base='G'),lineage2,comas09),
    Test(SNP(genome=ancestor,pos= 301341,base='A'),lineage3,stucki12),
    Test(SNP(genome=ancestor,pos=4266647,base='G'),lineage3,stucki12),
    Test(SNP(genome=ancestor,pos= 157129,base='T'),lineage3,comas09),
    Test(SNP(genome=ancestor,pos=3326554,base='A'),lineage4,stucki12), # same as in ancestor
    Test(SNP(genome=ancestor,pos=2154724,base='C'),lineage4,stucki12), # same as in ancestor
    Test(SNP(genome=ancestor,pos= 648856,base='T'),lineage4,stucki12), # same as in ancestor
    Test(SNP(genome=ancestor,pos=1377185,base='G'),lineage5,stucki12),
    Test(SNP(genome=ancestor,pos= 801959,base='T'),lineage5,stucki12),
    Test(SNP(genome=ancestor,pos=2859147,base='T'),lineage5,stucki12),
    Test(SNP(genome=ancestor,pos=2427828,base='C'),lineage6,stucki12),
    Test(SNP(genome=ancestor,pos= 378404,base='A'),lineage6,stucki12),
    Test(SNP(genome=ancestor,pos=4269522,base='A'),lineage6,stucki12),
    Test(SNP(genome=ancestor,pos=  14806,base='C'),lineage7,stucki12),
    Test(SNP(genome=ancestor,pos=1663221,base='G'),lineage7,stucki12),
    Test(SNP(genome=ancestor,pos= 497126,base='A'),lineage7,stucki12),
    Test(SNP(genome=ancestor,pos= 10819,base='T'),lineage9,zwyer21),
    Test(SNP(genome=ancestor,pos= 1810476,base='A'),lineage9,zwyer21),
    Test(SNP(genome=ancestor,pos= 24409,base='A'),lineage9,zwyer21),
    Test(SNP(genome=ancestor,pos= 3626025,base='A'),lineage9,zwyer21),
    Test(SNP(genome=ancestor,pos= 623728,base='T'),lineage9,zwyer21),
Test(SNP(genome=ancestor,pos= 17333,base='C'),lineage8,ngabonziza20),
    Test(SNP(genome=ancestor,pos= 103128,base='G'),lineage8,ngabonziza20),
    Test(SNP(genome=ancestor,pos= 162637,base='C'),lineage8,ngabonziza20),
    Test(SNP(genome=ancestor,pos= 2505582,base='T'),lineage8,ngabonziza20),
    Test(SNP(genome=ancestor,pos= 442577,base='C'),lineage8,ngabonziza20)]

sublineage_SNPs = [Test(SNP(genome=ancestor,pos=118468,base='T'),unk4_unk5_eu2,zwyer21),
Test(SNP(genome=ancestor,pos=2963597,base='A'),unk4_unk5_eu2,zwyer21),
Test(SNP(genome=ancestor,pos=540940,base='C'),unk4_unk5_eu2,zwyer21),
Test(SNP(genome=ancestor,pos=772740,base='T'),unk4_unk5_eu2,zwyer21),
Test(SNP(genome=ancestor,pos=2715125,base='C'),eu1_unk6_unk7,zwyer21),
Test(SNP(genome=ancestor,pos=288952,base='T'),eu1_unk6_unk7,zwyer21),
Test(SNP(genome=ancestor,pos=10727,base='G'),eu1_unk6_unk7,zwyer21),
Test(SNP(genome=ancestor,pos=2853155,base='A'),eu1_unk6_unk7,zwyer21),
Test(SNP(genome=ancestor,pos=3033644,base='C'),unk2,zwyer21),
Test(SNP(genome=ancestor,pos=1985034,base='C'),unk2,zwyer21),
Test(SNP(genome=ancestor,pos=4309281,base='T'),unk2,zwyer21),
Test(SNP(genome=ancestor,pos=318068,base='T'),unk2,zwyer21),
Test(SNP(genome=ancestor,pos=23297,base='T'),af2,zwyer21),
Test(SNP(genome=ancestor,pos=67076,base='T'),af2,zwyer21),
Test(SNP(genome=ancestor,pos=4193069,base='A'),af2,zwyer21),
Test(SNP(genome=ancestor,pos=3244324,base='G'),af2,zwyer21),
Test(SNP(genome=ancestor,pos=438444,base='C'),af1,zwyer21),
Test(SNP(genome=ancestor,pos=2848135,base='A'),af1,zwyer21),
Test(SNP(genome=ancestor,pos=2655824,base='C'),af1,zwyer21),
Test(SNP(genome=ancestor,pos=3072872,base='A'),af1,zwyer21),
Test(SNP(genome=ancestor,pos=2339255,base='A'),PZA_sus,zwyer21),
Test(SNP(genome=ancestor,pos=4155799,base='A'),PZA_sus,zwyer21),
Test(SNP(genome=ancestor,pos=2271555,base='A'),PZA_sus,zwyer21),
Test(SNP(genome=ancestor,pos=3581187,base='A'),PZA_sus,zwyer21),
Test(SNP(genome=ancestor,pos=4357503,base='G'),PZA_sus,zwyer21),
Test(SNP(genome=ancestor,pos=1065611,base='A'),unk3,zwyer21),
Test(SNP(genome=ancestor,pos=2428078,base='G'),unk3,zwyer21),
Test(SNP(genome=ancestor,pos=471576,base='C'),unk3,zwyer21),
Test(SNP(genome=ancestor,pos=2183681,base='A'),unk3,zwyer21),
Test(SNP(genome=ancestor,pos=3849394,base='A'),unk3,zwyer21),
Test(SNP(genome=ancestor,pos=2669115,base='A'),unk9,zwyer21),
Test(SNP(genome=ancestor,pos=594057,base='A'),unk9,zwyer21),
Test(SNP(genome=ancestor,pos=51459,base='C'),unk9,zwyer21),
Test(SNP(genome=ancestor,pos=3038566,base='A'),unk9,zwyer21),
Test(SNP(genome=ancestor,pos=1552837,base='T'),unk9,zwyer21)]

subsublineage_SNPs = [Test(SNP(genome=ancestor,pos=25306,base='G'),eu1,zwyer21),
Test(SNP(genome=ancestor,pos=2890668,base='T'),eu1,zwyer21),
Test(SNP(genome=ancestor,pos=4030259,base='G'),eu1,zwyer21),
Test(SNP(genome=ancestor,pos=3347070,base='T'),eu1,zwyer21),
Test(SNP(genome=ancestor,pos=1534929,base='T'),eu1,zwyer21),
Test(SNP(genome=ancestor,pos=1765767,base='C'),unk4,zwyer21),
Test(SNP(genome=ancestor,pos=2807764,base='A'),unk4,zwyer21),
Test(SNP(genome=ancestor,pos=1401050,base='T'),unk4,zwyer21),
Test(SNP(genome=ancestor,pos=1158434,base='A'),unk4,zwyer21),
Test(SNP(genome=ancestor,pos=636624,base='T'),unk4,zwyer21),
Test(SNP(genome=ancestor,pos=3720207,base='T'),eu2,zwyer21),
Test(SNP(genome=ancestor,pos=313605,base='A'),eu2,zwyer21),
Test(SNP(genome=ancestor,pos=2553395,base='A'),eu2,zwyer21),
Test(SNP(genome=ancestor,pos=4114819,base='C'),eu2,zwyer21),
Test(SNP(genome=ancestor,pos=2702178,base='A'),eu2,zwyer21),
Test(SNP(genome=ancestor,pos=2989191,base='C'),unk7,zwyer21),
Test(SNP(genome=ancestor,pos=2594632,base='A'),unk7,zwyer21),
Test(SNP(genome=ancestor,pos=4297210,base='G'),unk7,zwyer21),
Test(SNP(genome=ancestor,pos=4204696,base='A'),unk7,zwyer21),
Test(SNP(genome=ancestor,pos=1358422,base='T'),unk7,zwyer21),
Test(SNP(genome=ancestor,pos=2948297,base='C'),unk5,zwyer21),
Test(SNP(genome=ancestor,pos=1978431,base='T'),unk5,zwyer21),
Test(SNP(genome=ancestor,pos=1020732,base='C'),unk5,zwyer21),
Test(SNP(genome=ancestor,pos=2291557,base='G'),unk5,zwyer21),
Test(SNP(genome=ancestor,pos=1427067,base='G'),unk5,zwyer21),
Test(SNP(genome=ancestor,pos=220104,base='T'),unk6,zwyer21),
Test(SNP(genome=ancestor,pos=1288126,base='C'),unk6,zwyer21),
Test(SNP(genome=ancestor,pos=2243435,base='A'),unk6,zwyer21),
Test(SNP(genome=ancestor,pos=3676093,base='T'),unk6,zwyer21),
Test(SNP(genome=ancestor,pos=545403,base='T'),unk6,zwyer21),
Test(SNP(genome=ancestor,pos=128264,base='G'),unk2_BCG,zwyer21),
Test(SNP(genome=ancestor,pos=3423939,base='A'),unk2_BCG,zwyer21), 
Test(SNP(genome=ancestor,pos=2834424,base='A'),unk2_BCG,zwyer21), 
Test(SNP(genome=ancestor,pos=3095342,base='C'),unk2_BCG,zwyer21), 
Test(SNP(genome=ancestor,pos=2027617,base='A'),unk2_BCG,zwyer21)]
phylo = PhyloTestsuite(mainlineages_SNPs + sublineage_SNPs + subsublineage_SNPs, VERSION )