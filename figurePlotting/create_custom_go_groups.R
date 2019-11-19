#create_custom_go_groups.R

# PROTEIN FOLDING ---------------------------------------------------------

#build up protein folding df
foldingGos = c('GO:0006986',
                          'GO:0030968',
                          'GO:0034620',
                          'GO:0035966')
foldingNames = c('response to unfolded protein',
                            'endoplasmic reticulum unfolded protein response',
                            'cellular response to unfolded protein',
                            'response to topologically incorrect protein')


# PROTEIN DEGREDATION -----------------------------------------------------

proteinDegGos = c('GO:0010498',
                  'GO:0043161',
                  'GO:0032436',
                  'GO:1901800')
proteinDegNames = c('proteasomal protein catabolic process',
                    'proteasome-mediated ubiquitin-dependent protein catabolic process',
                    'positive regulation of proteasomal ubiquitin-dependent protein catabolic process',
                    'positive regulation of proteasomal protein catabolic process')


# CELL DEATH --------------------------------------------------------------

deathGos = c('GO:0010942',
             'GO:0012502',
             'GO:0012501')
deathNames = c('positive regulation of cell death',
               'induction of programmed cell death',
               'programmed cell death')

# RESPONSE TO OXIDATIVE STRESS --------------------------------------------

#build custom ros 
rosGos = c('GO:2000379',
           'GO:2000377',
           'GO:0042542',
           'GO:0000302')

rosNames = c('positive regulation of reactive oxygen species metabolic process',
             'regulation of reactive oxygen species metabolic process',
             'response to hydrogen peroxide',
             'response to reactive oxygen species')


# CELL CYCLE ARREST -------------------------------------------------------

arrestGos = c('GO:0071158',
              'GO:0044819',
              'GO:0072331',
              'GO:0030330')
arrestNames = c('positive regulation of cell cycle arrest',
                'mitotic G1/S transition checkpoint',
                'signal transduction by p53 class mediator',
                'DNA damage response, signal transduction by p53 class mediator')



# VESICLES ----------------------------------------------------------------

vesicleGos = c('GO:0006901',
               'GO:0006900',
               'GO:0006903')
vesicleNames = c('vesicle coating',
                 'vesicle budding from membrane',
                 'vesicle targeting')


# NF-kappaBeta ------------------------------------------------------------

nfkbGos = c('GO:0038061',
            'GO:0007250',
            'GO:0043122')
nfkbNames = c('NIK/NF-kappaB signaling',
              'activation of NF-kappaB-inducing kinase activity',
              'regulation of I-kappaB kinase/NF-kappaB signaling')


# ACTIVATION OF IMMUNE RESPONSE -------------------------------------------

immuneGos = c('GO:0002253')
immuneNames = c('activation of immune response')

# DNA REPLICATION ---------------------------------------------------------

dnaRepGos = c('GO:0006260',
              'GO:0030261',
              'GO:0071897')
dnaRepNames = c('DNA replication',
                'DNA biosynthetic process',
                'chromosome condensation')


# CELL DIVISION -----------------------------------------------------------

prolifGos = c('GO:0051785',
              'GO:0045840',
              'GO:0006312')
prolifNames = c('positive regulation of nuclear division',
                'positive regulation of mitotic nuclear division',
                'mitotic recombination')


# RIBOSOMES ---------------------------------------------------------------

#make custom for BP
ribosomeGos = c('GO:0042254')
ribosomeNames = c('ribosome biogenesis')



# CREATE CUSTOM DATAFRAME FOR OVERLAY -------------------------------------

upterms = c(foldingGos,
            proteinDegGos,
            deathGos,
            rosGos,
            arrestGos,
            vesicleGos,
            nfkbGos,
            immuneGos)

upnames = c(foldingNames,
            proteinDegNames,
            deathNames,
            rosNames,
            arrestNames,
            vesicleNames,
            nfkbNames,
            immuneNames)

upsummaries = c(rep('protein folding', times=length(foldingGos)),
                rep('protein degredation', times=length(proteinDegGos)),
                rep('cell death', times=length(deathGos)),
                rep('ROS', times=length(rosGos)),
                rep('cycle arrest', times=length(arrestGos)),
                rep('membrane vesicle', times=length(vesicleGos)),
                rep('NF-KB', times=length(nfkbGos)),
                rep('immune activition', times=length(immuneGos)))

upSelect = data.frame(term = upterms,
                    name = upnames,
                    summary = upsummaries,
                    stringsAsFactors=FALSE)


downterms = c(ribosomeGos,
              dnaRepGos,
              prolifGos)
downnames = c(ribosomeNames,
              dnaRepNames,
              prolifNames)
downsummaries = c(rep('ribosomes', times=length(ribosomeGos)),
                  rep('DNA replication', times=length(dnaRepGos)),
                  rep('cell division', times=length(prolifGos)))
downSelect = data.frame(term = downterms,
                        name = downnames,
                        summary = downsummaries,
                        stringsAsFactors = FALSE)



save(upSelect, downSelect, file='figurePlotting/selectGoGroups.Rdata')


#build a supplementary table to go with figure
sdat = rbind(upSelect, downSelect) %>% 
  select(summary, term, name)

sdat %>% 
  write_tsv(path='figurePlotting/customGoSets.tsv')







