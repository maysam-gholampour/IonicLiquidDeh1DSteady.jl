using IonicLiquidDeh1DSteady
using Documenter

DocMeta.setdocmeta!(IonicLiquidDeh1DSteady, :DocTestSetup, :(using IonicLiquidDeh1DSteady); recursive=true)

makedocs(;
    modules=[IonicLiquidDeh1DSteady],
    authors="maysam-gholampour <meysam.gholampoor@gmail.com> and contributors",
    sitename="IonicLiquidDeh1DSteady.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
