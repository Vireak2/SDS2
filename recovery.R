# Recovery
name = "5K_disassortative"
N_misclass = 8
sample_network = readRDS(paste("networks/",".rds", sep = name))
memberships =read.csv(paste("memberships/", ".csv", sep = name))
line_graph = line.graph(sample_network$graph)

membership_linkcomm = memberships$membership_linkcomm
membership_linkcommmod = memberships$membership_linkcommmod
membership_c = memberships$membership_c
membership_d = memberships$membership_d
membership_e = memberships$membership_e
membershipLS = memberships$membershipLS



mixed_membership_linkcomm = labelize(sample_network$node_edge, clean_membership(membership_linkcomm), mode = 'multi')
mixed_membership_linkcommmod = labelize(sample_network$node_edge, clean_membership(membership_linkcommmod), mode = 'multi')
mixed_membership_c = labelize(sample_network$node_edge, clean_membership(membership_c), mode='multi')
mixed_membership_d = labelize(sample_network$node_edge, clean_membership(membership_d), mode='multi')
mixed_membership_e = labelize(sample_network$node_edge, clean_membership(membership_e), mode='multi')
linkscan_membership = labelize(sample_network$node_edge, clean_membership(membershipLS), mode='multi')


colored_linkcomm = Color_Adjacency(sample_network$adjacency,
                                   sample_network$edges,
                                   clean_membership(membership_linkcomm))
colored_linkcommmod = Color_Adjacency(sample_network$adjacency,
                                      sample_network$edges,
                                      clean_membership(membership_linkcommmod))
colored_c = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            clean_membership(membership_c))
colored_d = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            clean_membership(membership_d))
colored_e = Color_Adjacency(sample_network$adjacency,
                            sample_network$edges,
                            clean_membership(membership_e))
colored_linkscan = Color_Adjacency(sample_network$adjacency,
                                   sample_network$edges,
                                   clean_membership(membershipLS))




# Modularity
mod_linkcomm = modularity(line_graph, membership_linkcomm+1)
mod_linkcommmod = modularity(line_graph, membership_linkcommmod+1)
mod_c = modularity(line_graph, membership_c)
mod_d = modularity(line_graph, membership_d)
mod_e = modularity(line_graph, membership_e)
mod_linkscan = modularity(line_graph, membershipLS)
modularities = c(mod_linkcomm,mod_linkcommmod, mod_c, mod_d, mod_e, mod_linkscan) 

# Expected
expected = expected_communities(sample_network, TRUE)
exp_linkcomm = norm(expected - number_communities(mixed_membership_linkcomm), type="2")
exp_linkcommmod = norm(expected - number_communities(mixed_membership_linkcommmod), type="2")
exp_c = norm(expected - number_communities(mixed_membership_c), type="2")
exp_d = norm(expected - number_communities(mixed_membership_d), type="2")
exp_e = norm(expected - number_communities(mixed_membership_e), type="2")
exp_linkscan = norm(expected - number_communities(linkscan_membership), type="2")
expectations = c(exp_linkcomm,exp_linkcommmod, exp_c, exp_d, exp_e, exp_linkscan) 

# misclassification approx
misclass_linkcomm = misclassification_edges(sample_network, colored_linkcomm, N_misclass)[[2]]
misclass_linkcommmod = misclassification_edges(sample_network, colored_linkcommmod, N_misclass)[[2]]
misclass_c = misclassification_edges(sample_network, colored_c, N_misclass)[[2]]
misclass_d = misclassification_edges(sample_network, colored_d, N_misclass)[[2]]
misclass_e = misclassification_edges(sample_network, colored_e, N_misclass)[[2]]
misclass_linkscan = misclassification_edges(sample_network, colored_linkscan, N_misclass)[[2]]
misclassifications = c(misclass_linkcomm,misclass_linkcommmod, misclass_c, misclass_d, misclass_e, misclass_linkscan)

# misclassifications on vertex using diagonal prediction
v_mis_LC = misclassification_vertex(sample_network, colored_linkcomm)
v_mis_LCmod = misclassification_vertex(sample_network, colored_linkcommmod)
v_mis_c = misclassification_vertex(sample_network, colored_c)
v_mis_d = misclassification_vertex(sample_network, colored_d)
v_mis_e = misclassification_vertex(sample_network, colored_e)
v_mis_LS = misclassification_vertex(sample_network, colored_linkscan)
v_missclassifications=c(v_mis_LC, v_mis_LCmod,v_mis_c, v_mis_d, v_mis_e, v_mis_LS) 

c(mod_linkscan, exp_linkscan, misclass_linkscan, v_mis_LS)

measurement = data.frame(Methods = c('Linkcomm','Linkcomm + mod', 'Imperial C','Imperial D', 'Imperial E', 'LinkScan*'),
                         Modularity = modularities,
                         Expectation_L2error = expectations,
                         Approx_misclassification = misclassifications,
                         vertex_misclassifications = v_missclassifications)
measurement

write.csv(measurement, file = paste("measures_csv/", ".csv", sep=name), row.names = FALSE)
