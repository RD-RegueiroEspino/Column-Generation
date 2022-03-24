using JuMP, LinearAlgebra, GLPK, DelimitedFiles

function init_A(N,L,l,n)
    A = zeros((N,N))
    for i in 1:N
        A[i,i]=min(floor(Int,L/l[i]),round(Int,n[i]))
    end
    return(A)
end

function print_matrix(f,M)
    nrow = size(M)[1]
    for i in 1:nrow
        println(f,M[i,:])
    end
end

function init_master(A,L,l,n,f,i)
    N = length(l)
    ncols = size(A)[2]
    model = Model(GLPK.Optimizer)
    @variable(model, lambda[1:ncols]>=0)
    @constraint(model, c[j=1:N], sum(A[j,p]*lambda[p] for p in 1:ncols)>=n[j])
    @objective(model, Min, sum(lambda[p] for p in 1:ncols))
    println(f, "RMP $i\n")
    println(f, model)
    optimize!(model)
    println(f, "The solution is:")
    res = zeros(ncols)
	for i in 1:ncols
        res[i] = JuMP.value(lambda[i])
		println(f, "lambda[", i, "] = ", string(JuMP.value(lambda[i])))
	end
    println(f,"\nThe value of the dual variables are:")
    dual_coeff = zeros(N)
    for i in 1:length(l)
        dual_coeff[i]=dual(c[i])
        println(f, "pi[", i, "] = ", dual_coeff[i])
    end
    println(f)
    return res,dual_coeff
end

function auxiliary_problem(mp,l,L,f,i)
    N = length(l)
    ncols = length(l)
    println(f, "AP $i\n")
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:ncols]>=0, Int)
    @constraint(model, sum(l[p]*x[p] for p in 1:ncols)<=L)
    @objective(model, Max, sum(mp[p]*x[p] for p in 1:ncols))
    println(f, model)
    optimize!(model)
    x_values = zeros(ncols)
    println(f, "The solution is:")
	for i in 1:ncols
        x_values[i]=value(x[i])
		println(f, "x[", i, "] = ", string(JuMP.value(x[i])))
	end
    println(f)
    return x_values
end

function resolution(A,L,l,n,f)
    i=1
    reduced_cost = -1
    new_column = 0
    res = 0
    while (reduced_cost<0)
        res, mp = init_master(A,L,l,n,f,i) 
        new_column = auxiliary_problem(mp,l,L,f,i) 
        reduced_cost = 1 - sum(new_column[p]*mp[p] for p in 1:length(new_column))
        println(f, "Hence, the reduced cost is $reduced_cost\n")
        if reduced_cost < 0
            A = hcat(A,new_column)
            println(f,"So we have a new column a$(size(A)[2])")
            println(f,new_column)
            println(f)
        end
        if reduced_cost >= 0 
            println(f,"Thus we do not have a new column
and the optimal solution was found at RMP $i\n")
        end
        i+=1
    end
    return res,A
end

function presentation_exercice(N,L,l,n,A,f)
    println(f, "A company produces steel bars with L = $L m and cuts the bars
for the costumers according to their necessities.
Now, the company has to satisfy the following demand:\n")
    println(f, "l_i | number of pieces needed")
    for i in 1:N
        println(f, "$(l[i])  | $(n[i])")
    end
    println(f)
    println(f, "The costs c_p are all equal to one.
We want to minimize the number of steel bars we need to cut
in order to satisfy the demand.\n")
    println(f,"In order to find an initial solution,
we compute how many pieces of each length fits in one bar.\n")
    println(f, "So our matrix A at the first iteration is:")
    print_matrix(f,A)
    println(f)
end

function conclusion_exercice(A,l,res,f)
    nrow = size(A)[1]
    for i in 1:length(res)
        res[i]=round(ceil(res[i]), digits=0)
    end
    println(f,"Rounding the solution of the ultimate RMP we have the following.\n")
    for i in 1:length(res)
        if res[i]>0
            println(f, "- cut $(res[i]) entire bar(s) in pattern $i, which means")
            for j in 1:(size(A)[1])
                if A[j,i]>0
                    println(f, "    $(A[j,i]) piece(s) of $(l[j]) m")
                end
            end
        end
    end
    println(f,"\nThis solution gives us :")
    prod = A*res
    for i in 1:nrow
        println(f,"- $(prod[i]) piece(s) of $(l[i])m")
    end
end

function ex1()
    L=10
    l=[2,3,4]
    n=[6,4,2]
    N=length(l)
    A = init_A(N,L,l,n)
    f = open("output_ex1.txt", "w")
    println(f,"EX1\n")
    presentation_exercice(N,L,l,n,A,f)
    res, A = resolution(A,L,l,n,f)
    conclusion_exercice(A,l,res,f)
    close(f)
end

ex1()

function ex2()
    L=100
    l=[22,42,52,53,78]
    n=[45,38,25,11,12]
    N=length(l)
    A = init_A(N,L,l,n)
    f = open("output_ex2.txt", "w")
    println(f,"EX2\n")
    presentation_exercice(N,L,l,n,A,f)
    res, A = resolution(A,L,l,n,f)
    conclusion_exercice(A,l,res,f)
    close(f)
end

ex2()