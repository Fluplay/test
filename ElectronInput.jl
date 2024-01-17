
let
    ####################
    # Global Variables #
    ####################

    N = -1
    J = 1

    ########################################################
    # construct site and Hamiltonian (for the Spin-1 case) #
    ########################################################

    sites = siteinds("S=1",N)


    os = OpSum()

    for j = 1 : N-1
        os += J,"Sx",j,"Sx",j+1
        os += J,"Sy",j,"Sy",j+1
        os += J,"Sz",j,"Sz",j+1
    end

    H = MPO(os,sites)
    println("")

    ###########################################
    # Hyper parameters for DMRG calculations  #
    ###########################################

    nsweeps = 10
    maxdim = [20,80,120,240,480,600,800,1000]
    mindim = 80
    cutoff = 4E-10

    ###########################################
    # DMRG Calculation with random state psi0 #
    ###########################################

    psi0 = randomMPS(sites)
    energy,psi = dmrg(H,psi0 ;nsweeps,mindim,maxdim,cutoff,outputlevel = 1)

    # using ITensors.HDF5
    # f = h5open("Heisenberg.h5","w")
    # write(f,"T",psi)
    # close(f)


    ########################################
    # Function for Trotter gate O(delta^3) #
    ########################################
    function MAKEGATE(_tau,_sites,J)
        gatesLRRL = ITensor[]
        N = length(_sites)
        for i = 1:N-1
            j = i+1
            si = _sites[i]
            sj = _sites[j]
            gatesLR = ITensor[]
            hj = J*op("Sx",si)*op("Sx",sj) + J*op("Sy",si)*op("Sy",sj) + J*op("Sz",si)*op("Sz",sj)
            Gj = exp(-im * _tau * (hj/2))
            push!(gatesLRRL,Gj)
        end

        for i = N-1:-1:1
            j = i+1
            si = _sites[i]
            sj = _sites[j]
            gatesLR = ITensor[]
            hj = J*op("Sx",si)*op("Sx",sj) + J*op("Sy",si)*op("Sy",sj) + J*op("Sz",si)*op("Sz",sj)
            Gj = exp(-im * _tau * (hj/2))
            push!(gatesLRRL,Gj)
        end
        gatesRLLR = reverse(gatesLRRL)
        return gatesLRRL,gatesRLLR
    end

    gatesLRRL,gatesRLLR =  MAKEGATE(0.1,sites,1)

    #################################
    # Make Mixture for Ground state #
    #################################

    ME = op("S+",sites[N])
    MF = op("S+",sites[1])
   #MC = op("Sz",sites[Int(N/2)])
    psi = apply(ME,psi)
    psi = apply(MF,psi)
   # psi = apply(MC,psi)


-- INSERT --                                        