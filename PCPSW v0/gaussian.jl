#computationalthinking.mit.edu
#Gaussian Distribution
#Probability Density Function(pdf)
#Cumulative Density Function(cdf)

using Statistics
#using Plots
using SpecialFunctions

begin
    abstract type RandomVariable end
    abstract type DiscreteRandomVariable <: RandomVariable end
    abstract type ContinousRandomVariable <: RandomVariable end
end

begin
    struct Gaussian <: ContinousRandomVariable
        μ   #mean
        σ   #standard deviation
    end
    Gaussian() = Gaussian(0.0, 1.0) #normalised Gaussian with mean 0.0 adn variance 1.0
end

begin
    Statistics.mean(X::Gaussian) = X.μ
    Statistics.std(X::Gaussian) = X.σ
end

#variance
Statistics.var(X::Gaussian) = std(X)^2

#Probability distribution of a Gaussian (Normal)
pdf_gaussian(X::Gaussian) = x-> exp(-0.5 * ((x - X.μ)^2 / X.σ^2)) / √(2π * X.σ^2)
cdf_gaussian(X::Gaussian) = x-> 0.5 * (erf((x - X.μ) / √(2 * X.σ^2)) + 1)

#Log Normal
pdf_LogNorm(X::Gaussian) = x-> exp(-0.5 * ((log(x) - X.μ)^2 / X.σ^2)) / (x * √(2π * X.σ^2))
cdf_LogNorm(X::Gaussian) = x-> 0.5 * (erf((log(x) - X.μ) / √(2 * X.σ^2)) + 1)

#example using Gaussian distribution
begin
#    plot(pdf_gaussian(Gaussian(0.0, 1.0)), leg=false)
#    xlims!(-6, 6)
#    ylims!(0, 0.5)
end

begin
#    g = Gaussian(4.08, 4.08 * 0.3125)
#    plot((1:75*12)/100, pdf_gaussian(g))
#    plot((1:75*12)/100, cdf_gaussian(g))
#    plot((1:75*12), cdf_LogNorm(g))
#    plot((1:75*12), pdf_LogNorm(g))
end

begin
#    plot(cdf_gaussian(Gaussian(0.0, 1.0)), leg=false)
#    xlims!(-6, 6)
#    ylims!(0, 1.2)
end

#Gaussian distribution
Base.rand(X::Gaussian) = X.μ + √(X.σ^2) * randn()

#example histogram Gaussian distribution
#histogram!([rand(Gaussian(0.0, 1.0)) for i in 1:10^4], alpha=0.5, norm=true)
