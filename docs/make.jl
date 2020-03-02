using Documenter, RANSAC

makedocs(sitename="RANSAC.jl",
        pages=["Home" => "index.md",
                "Example" => "example.md",
                "Efficient RANSAC" => "ransac.md",
                "API" => "api.md"])


deploydocs(repo="github.com/cserteGT3/RANSAC.jl.git")
