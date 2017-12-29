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
    parameter1 = set_parameter_defaults(keyword1)

    @test parameter1.travel == "grid"
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

    # first test "genetic counseling 1" files
    par = parameter.par #this is [0.0]
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
        person, nuclear_family, locus, keyword)

    # note the following probabilities are provided in the locus frame
    dic_locus1 = Dict("1" => 0.62, "2" => 0.17, "3" => 0.14, "4" => 0.07)
    dic_locus2 = Dict("1" => 0.003, "2" => 0.997)
    locus_3_sum = 0.41 + 0.14 + 0.03 + 0.39 #need to normalize allele freq to 1 at locus 3
    dic_locus3 = Dict("1" => 0.41/locus_3_sum, "2" => 0.14/locus_3_sum,
        "3" => 0.03/locus_3_sum, "4" => 0.39/locus_3_sum)

    for n = instruction.start[1]:instruction.finish[1]-1
        # n is the variable iterating from instruction.start[ped] to instruction.finish[ped].
        operation = instruction.operation[n]
        start = instruction.extra[n][1]
        finish = instruction.extra[n][2]
        i = instruction.extra[n][3]
        if operation != penetrance_and_prior_array continue end #this avoids some array access errors
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

@testset "penetrance function" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["gender_neutral"] = true
    keyword["lod_score_table"] = "Lod_Score_Frame.txt"
    # keyword["eliminate_genotypes"] = true
    # keyword["lump_alleles"] = true
    keyword["parameters"] = 1
    keyword["points"] = 9
    process_keywords!(keyword, "two-point linkage Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
        locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    parameter = set_parameter_defaults(keyword)
end

@testset "transmission function" begin
    
end

@testset "wrapper and basics" begin

end


