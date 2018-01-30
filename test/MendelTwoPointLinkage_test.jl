using MendelTwoPointLinkage
using MendelBase
using Search
using SearchSetup
using DataFrames

@testset "initialize_optimization_two_point_linkage" begin
    keyword1 = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword1["gender_neutral"] = true
    keyword1["lod_score_table"] = "Lod_Score_Frame.txt"
    keyword1["parameters"] = 1
    keyword1["points"] = 9
    process_keywords!(keyword1, "two-point linkage Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
        locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword1)
    parameter1 = set_parameter_defaults(keyword1) #must process_keywords before setting parameters

    # not yet called initialize_optimization
    @test parameter1.travel == "grid" #from control file
    @test size(parameter1.grid) == (9, 1) #(points, par)
    @test all(parameter1.grid .== 0.0)
    @test size(parameter1.par) == (1,)
    @test parameter1.par[1] == 0.0
    @test size(parameter1.name) == (1,)
    @test size(parameter1.min) == (1,)
    @test size(parameter1.max) == (1,)
    @test parameter1.min[1] == -Inf
    @test parameter1.max[1] == Inf
    @test parameter1.parameters == 1

    parameter1 = MendelTwoPointLinkage.initialize_optimization_two_point_linkage!(
        locus, parameter1, keyword1)

    @test parameter1.travel == "grid"
    @test size(parameter1.grid) == (9, 1) #(points, par)
    @test all(parameter1.grid .== [0.5; 0.4; 0.3; 0.2; 0.15; 0.1; 0.05; 0.01; 0.001])
    @test size(parameter1.par) == (1,)
    @test parameter1.par[1] == 0.0
    @test size(parameter1.name) == (1,)
    @test size(parameter1.min) == (1,)
    @test size(parameter1.max) == (1,)
    @test parameter1.min[1] == -Inf
    @test parameter1.max[1] == Inf
    @test parameter1.parameters == 1


    keyword2 = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword2["gender_neutral"] = true
    keyword2["lod_score_table"] = "Lod_Score_Frame.txt"
    keyword2["parameters"] = 2
    keyword2["points"] = 9
    process_keywords!(keyword2, "two-point linkage Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
        locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword2)
    keyword2["travel"] = "search"

    parameter2 = set_parameter_defaults(keyword2)
    parameter2 = MendelTwoPointLinkage.initialize_optimization_two_point_linkage!(
        locus, parameter2, keyword2)

    @test parameter2.travel == "search"
    @test size(parameter2.par) == (2,)
    @test all(parameter2.par .== 0.5 - 1e-5)
    @test size(parameter2.name) == (2,)
    @test parameter2.name[1] == "xxtheta"
    @test parameter2.name[2] == "xytheta"    
    @test size(parameter2.min) == (2,)
    @test size(parameter2.max) == (2,)
    @test all(parameter2.min .== 1e-5)
    @test all(parameter2.max .== 0.5)
    @test parameter2.parameters == 2
end

@testset "prior function" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gender_neutral"] = true
    keyword["lod_score_table"] = "Lod_Score_Frame.txt"
    keyword["parameters"] = 1
    keyword["points"] = 9
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "two-point linkage Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    keyword["eliminate_genotypes"] = true
    keyword["lump_alleles"] = true

    #
    # first test "genetic counseling 1" files
    #
    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    # use dictionaries to assign prior probabilities
    # the numbers are provided in the locus frame
    dic_locus1 = Dict("1" => 0.62, "2" => 0.17, "3" => 0.14, "4" => 0.07)
    dic_locus2 = Dict("1" => 0.003, "2" => 0.997)
    locus_3_sum = 0.41 + 0.14 + 0.03 + 0.39 #need to normalize allele freq to 1 at locus 3
    dic_locus3 = Dict("1" => 0.41/locus_3_sum, "2" => 0.14/locus_3_sum,
        "3" => 0.03/locus_3_sum, "4" => 0.39/locus_3_sum)

    # loop through to compute probabilities
    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    for ped in 1:pedigree.pedigrees
        for n = instruction.start[ped]:instruction.finish[ped]-1
            operation = instruction.operation[n]
            start = instruction.extra[n][1]
            finish = instruction.extra[n][2]
            i = instruction.extra[n][3]
            if operation != penetrance_and_prior_array continue end #avoids some array access errors
            if person.mother[i] != 0 continue end #prior prob doesnt exist for non founder

            #
            # Construct the parent's multiple locus genotypes.
            #
            genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
            multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                                genotypes, i)

            for j = 1:genotypes
                prob = MendelTwoPointLinkage.prior_two_point_linkage(person, locus, 
                    multi_genotype[:, :, j], par, keyword, start, finish, i)
                answer = 1.0

                # tally probabilities contributed by the 3 locus 
                answer *= dic_locus1[string(multi_genotype[1, 1, j])]
                answer *= dic_locus1[string(multi_genotype[2, 1, j])]

                answer *= dic_locus2[string(multi_genotype[1, 2, j])]
                answer *= dic_locus2[string(multi_genotype[2, 2, j])]

                answer *= dic_locus3[string(multi_genotype[1, 3, j])]
                answer *= dic_locus3[string(multi_genotype[2, 3, j])]

                @test answer == prob
            end
        end
    end
end

@testset "penetrance function" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gender_neutral"] = true
    keyword["lod_score_table"] = "Lod_Score_Frame.txt"
    keyword["parameters"] = 1
    keyword["points"] = 9
    process_keywords!(keyword, "two-point linkage Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
        locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    parameter = set_parameter_defaults(keyword)

    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)
    par = parameter.par #this is [0.0]

    # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
    for ped in 1:pedigree.pedigrees
        for n = instruction.start[ped]:instruction.finish[ped]-1
            operation = instruction.operation[n]
            start = instruction.extra[n][1]
            finish = instruction.extra[n][2]
            i = instruction.extra[n][3]
            if operation != penetrance_and_prior_array continue end #this avoids some array access errors

            #
            # Construct the parent's multiple locus genotypes.
            #
            genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
            multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
                                                genotypes, i)

            for j = 1:genotypes
                number = MendelTwoPointLinkage.penetrance_two_point_linkage(person, locus, 
                    multi_genotype[:, :, j], par, keyword, start, finish, i)
                @test number == 1.0 #this is always 1 
            end
        end
    end
end

@testset "transmission function" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gender_neutral"] = true
    keyword["lod_score_table"] = "Lod_Score_Frame.txt"
    keyword["parameters"] = 1
    keyword["points"] = 9
    parameter = set_parameter_defaults(keyword)
    process_keywords!(keyword, "two-point linkage Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    for ped in 1:pedigree.pedigrees
	    for n = instruction.start[ped]:instruction.finish[ped]-1
	        operation = instruction.operation[n]
	        if operation != transmission_array; continue; end #this avoids some array access errors

	        start = instruction.extra[n][1]
	        finish = instruction.extra[n][2]
	        i = instruction.extra[n][3]
	        j = instruction.extra[n][4]

	        # need 2 genotypes to run transmission, so we construct them
	        # and give them the same names as used in elston_stewart_evaluation
	        i_genotypes = MendelBase.genotype_count(person, locus, i, start, finish)
	        multi_genotype = MendelBase.construct_multigenotypes(person, locus, start, finish,
	                                            i_genotypes, i)

	        maternal = !person.male[i]
	        j_genotypes = MendelBase.genotype_count(person, locus, j, start, finish)
	        gamete = MendelBase.construct_gametes(person, locus, start, finish, j_genotypes, j,
	                                     maternal)

	        # loop through all possible genotypes, and see if transmission probability is correct
	        for l = 1:j_genotypes
	            for k = 1:i_genotypes
	                trans = MendelTwoPointLinkage.transmission_two_point_linkage(person, 
	                    locus, gamete[:, l], multi_genotype[:, :, k], par, keyword, 
	                    start, finish, i, j)
	                #
	                #compute the probaiblity that a multi_genotype produce a certain gamete
	        		#0.35 is recombination fraction between locus 2 and 3 given by calling locus.theta
	        		#Note between locus 1 and 2 the recomb fraction is 0 since they are 0 morgans apart
	        		#
	        		recomb_12 = 0.0
		            recomb_23 = 0.35

	            	if !(gamete[1, l] in multi_genotype[:, 1, k]) ||
	            	   !(gamete[2, l] in multi_genotype[:, 2, k]) ||
	            	   !(gamete[3, l] in multi_genotype[:, 3, k])
	            	    #case1: an allele came out of nowhere
	            		@test trans == 0.0
	            	elseif all(multi_genotype[1, :, k] .== multi_genotype[2, :, k])
	            		#case2: if two multi-locus-genotype is identical, then gamete will be too
	            		@test trans == 1.0
	            	elseif all(multi_genotype[1, 1:2, k] .== multi_genotype[2, 1:2, k])
            			#case3: The first two allele of both multi locus genotype is the same.
            			#		Then the 3rd (which must be different) passes down with prob 0.5
						@test trans == 0.5
	       			elseif !(all(multi_genotype[1, 1:2, k] .== multi_genotype[2, 1:2, k])) 
						# When the first two allele of the two multi locus genotype is NOT the same, 
						# we need to break it up into different cases
						if !(all(gamete[1:2, l] .== multi_genotype[1, 1:2, k])) &&
						   !(all(gamete[1:2, l] .== multi_genotype[2, 1:2, k])) 
						    # case4: the first two alleles in gamete must be the same as one of the first two 
						    # alleles in the multi locus genotype, since recombination between locus 1 
						    # and 2 is zero. So the prob that they're not the same is 0 
						    @test trans == recomb_12 * 0.5
						elseif multi_genotype[1, 3, k] == multi_genotype[2, 3, k]
							# case5: if the 3rd locus in both multi-locus-genotype is the same, then one of the 
							# first two locus pass down with prob 0.5, and we sum the case where 
							# recombination happens or not since both will create the desired gamete. 
	       					@test trans == 0.5 * recomb_23 + 0.5 * (1 - recomb_23) 
	       				elseif all(gamete[:, l] .== multi_genotype[1, :, k]) ||
	       					   all(gamete[:, l] .== multi_genotype[2, :, k])
	       					# case6: the first two alleles from one of the multi-locus-genotype was passed
	       					# down, and recombination didn't happen so that the third locus is too
	       					@test trans == 0.5 * (1 - recomb_23) 
	       				else
	       					# case7: the first two alleles from one of the multi-locus-genotype was passed
	       					# down, and recombination did happen so that the third locus from the other 
	       					# one is passed down.
	       					@test trans == 0.5 * recomb_23
						end
	           		else
		                @test_throws(MethodError, "all cases should have been considered")
	       			end
	            end
	        end
	    end    
	end
end

@testset "wrapper and basics" begin
    #
    # There should be no error when running two_point_linkage_option and TwoPointLinkage, given the input. 
    # If there were, then program should halt and would not return "false" or "nothing"
    #
	keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gender_neutral"] = true
    keyword["lod_score_table"] = "Lod_Score_Frame.txt"
    keyword["parameters"] = 1
    keyword["points"] = 9
    process_keywords!(keyword, "two-point linkage Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    
    execution_error = MendelTwoPointLinkage.two_point_linkage_option(pedigree, person, nuclear_family,
      locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
    final_test = TwoPointLinkage("two-point linkage Control.txt")

    @test execution_error == false
    @test final_test == nothing
end


#coverage = ?









