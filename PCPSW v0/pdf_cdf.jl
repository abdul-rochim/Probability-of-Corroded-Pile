using Distributions
using Statistics

pdfNorm(mu, sig, x) = @. pdf(Normal(mu,sig),x)
pdfLogNorm(mu, sig, x) = @. pdf(LogNormal(mu,sig),x)

cdfNorm(mu, sig, x) = @. cdf(Normal(mu,sig),x)
cdfLogNorm(mu, sig, x) = @. cdf(LogNormal(mu,sig),x)

@unpack pdf_which_Ccr, pdf_which_C = Params_Data()
if pdf_which_Ccr == 1
    pdf_used_Ccr = pdfNorm
elseif pdf_which_Ccr == 2
    pdf_used_Ccr = pdfLogNorm
else
    println("Choose which pdf(Ccr) you want to use!")
end

if pdf_which_C == 1
    pdf_used_C = pdfNorm
elseif pdf_which_C == 2
    pdf_used_C = pdfLogNorm
else
    println("Choose which pdf(C) you want to use!")
end

#example
begin
    #plot((1:75*12)/100, pdfNorm(4.08, 4.08*0.3125, (1:75*12)/100))
    #plot((1:75*12), pdfLogNorm(4.08, 4.08*0.3125, (1:75*12)))
end