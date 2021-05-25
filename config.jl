# config.jl for RCCP
# Prepares environment to run:

#using JuMP, MathProgBase, Gurobi
solver = "CPLEX"
if solver == "Gurobi"
  using JuMP, Gurobi
elseif solver == "GLPKLP" || solver == "GLPKMIP"
  using JuMP, GLPKMathProgInterface
elseif solver == "CPLEX"
  using  JuMP, CPLEX  #,JuMPeR
else
  println("NO SOLVER DEFINED")
end

MACHINE = "laptop"
if MACHINE == "laptop"
  if Sys.isapple()
    home_prefix = "/Users/levorato"
    git_prefix = "doutorado"
  else
    home_prefix = "/projetos/czt0"
    git_prefix = "doutorado.old"
  end
else
  home_prefix = joinpath("C:\\", "Users", "Public")
  git_prefix = "doutorado"
end
antoine_instances_folder = joinpath(home_prefix, git_prefix, "robusto", "RCCP", "instances", "antoine_skew")
toy_instances_folder = joinpath(home_prefix, git_prefix, "robusto", "RCCP", "instances", "toy")
instances_folder = joinpath(home_prefix, git_prefix,  "robusto", "RCCP", "instances")
project_folder = joinpath(home_prefix, git_prefix, "robusto", "RCCP")
output_folder = joinpath(home_prefix, "rccp_experiments")
println("*** Project folder is $(project_folder)")
println("*** Instances folder is $(instances_folder)")
println("*** Output folder is $(output_folder)")

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
