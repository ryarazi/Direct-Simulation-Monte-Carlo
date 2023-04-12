using Plots
using Printf
using Distributions
using HypothesisTests


function animate_system(samp)
    anim = @animate for i in 1:length(samp.t)
        sys = sys[i]
        scatter(
            sys.r[1, :], sys.r[2, :], #sys.r[3, :],
            xlim = (0, sys.L[1]), ylim = (0, sys.L[2]), zlim = (0, sys.L[3]),
            legend = false,
            marker = 3,
        )
    end
    gif(anim, "anim.gif", fps = 5)
end

function plot_thermalization(samp)
    Nsim = samp.sys[1].Nsim

    E_total = sum(samp.sysProp[end].E)
    kT = (2/3) * (E_total/Nsim) # because N*E_mean = 3/2 * kT
    E_array = range(0.05*kT, 6*kT, 200)
    
    MB = Gamma(3/2, kT)
    ymax = maximum(pdf.(MB, E_array))*1.3

    anim = @animate for i in 1:length(samp.t)
        AD_pvalue = pvalue(OneSampleADTest(samp.sysProp[i].E, MB))
        title = @sprintf "time=%.2g, Anderson-Darling p value=%.2f" samp.t[i] AD_pvalue
        histogram(samp.sysProp[i].E, bins=range(0.05*kT, 6*kT, 30), title=title, ylim=(0, ymax), normalize=:pdf, label="Particle Energy")
        plot!(E_array, pdf.(MB, E_array), xlim=(0.05*kT, 6*kT), ylim=(0, ymax),label="Maxwell-Boltzmann")
    end

    gif(anim, fps = 5)
end