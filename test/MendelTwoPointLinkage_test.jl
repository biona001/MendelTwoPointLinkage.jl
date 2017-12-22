using MendelTwoPointLinkage
using MendelBase
using Search
using SearchSetup
using DataFrames  

# function return_keyword()
#     keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
#     keyword["gender_neutral"] = true
#     keyword["lod_score_table"] = "Lod_Score_Frame.txt"
#     keyword["parameters"] = 1
#     keyword["points"] = 9
#     process_keywords!(keyword, "two-point linkage Control.txt", "")
#     (pedigree, person, nuclear_family, locus, snpdata,
#         locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
#         read_external_data_files(keyword)
#     return keyword
# end

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

@testset "penetrance function" begin
    
end

@testset "prior function" begin
    
end

@testset "transmission function" begin
    
end

@testset "wrapper and basics" begin

end


