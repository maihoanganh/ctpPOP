
function test_OPF_problem_case3_lmbd()
    data = PowerModels.parse_file("../test/opf_test/pglib_opf_case3_lmbd.m")
    
    run_OPF(data)

    println("====================")
    
end

function test_OPF_problem_case89_pegase__api()
    data = PowerModels.parse_file("../test/opf_test/pglib_opf_case89_pegase__api.m")
    
    run_OPF(data)

    println("====================")
    
end
