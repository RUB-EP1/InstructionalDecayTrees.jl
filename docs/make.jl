using InstructionalDecayTrees
using Documenter

const DOCS = @__DIR__
const QUARTO_PAGES = [
    (
        qmd = joinpath(DOCS, "forked_execution.qmd"),
        gfm = joinpath(DOCS, "forked_execution.md"),
        md = joinpath(DOCS, "src", "forked_execution.md"),
        id = "forked-execution",
        edit = "../forked_execution.qmd",
    ),
    (
        qmd = joinpath(DOCS, "wigner_su2_so3.qmd"),
        gfm = joinpath(DOCS, "wigner_su2_so3.md"),
        md = joinpath(DOCS, "src", "wigner_tutorial.md"),
        id = "wigner",
        edit = "../wigner_su2_so3.qmd",
    ),
    (
        qmd = joinpath(DOCS, "massless_wigner_limit.qmd"),
        gfm = joinpath(DOCS, "massless_wigner_limit.md"),
        md = joinpath(DOCS, "src", "massless_wigner_limit.md"),
        id = "massless-wigner-limit",
        edit = "../massless_wigner_limit.qmd",
    ),
    (
        qmd = joinpath(DOCS, "helicity_frame_conventions.qmd"),
        gfm = joinpath(DOCS, "helicity_frame_conventions.md"),
        md = joinpath(DOCS, "src", "helicity_frame_conventions.md"),
        id = "helicity-frame-conventions",
        edit = "../helicity_frame_conventions.qmd",
    ),
]

function render_quarto_page!(page)
    cd(DOCS) do
        run(`quarto render $(basename(page.qmd)) --to gfm`)
    end
    isfile(page.gfm) || error("expected Quarto output at $(page.gfm)")
end

function documenter_tutorial_page(page)
    gfm_path = page.gfm
    body = read(gfm_path, String)
    body = replace(
        body,
        r"^# ([^\r\n]+)"m => SubstitutionString("# [\\1](@id $(page.id))");
        count = 1,
    )
    meta = "```@meta\nCurrentModule = InstructionalDecayTrees\nEditURL = \"$(page.edit)\"\n```\n\n"
    return meta * body
end

DocMeta.setdocmeta!(
    InstructionalDecayTrees,
    :DocTestSetup,
    :(using InstructionalDecayTrees);
    recursive = true,
)

for page in QUARTO_PAGES
    render_quarto_page!(page)
    write(page.md, documenter_tutorial_page(page))
end

makedocs(;
    modules = [InstructionalDecayTrees],
    authors = "Mikhail Mikhasenko and contributors",
    repo = "https://github.com/RUB-EP1/InstructionalDecayTrees.jl/blob/{commit}{path}#{line}",
    sitename = "InstructionalDecayTrees.jl",
    doctest = true,
    checkdocs = :exports,
    format = Documenter.HTML(;
        canonical = "https://rub-ep1.github.io/InstructionalDecayTrees.jl",
        repolink = "https://github.com/RUB-EP1/InstructionalDecayTrees.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Forked execution" => "forked_execution.md",
        "API reference" => "api.md",
        "Wigner angles: SO(3) vs SU(2)" => "wigner_tutorial.md",
        "Particle-2 frame conventions" => "helicity_frame_conventions.md",
        "Small-mass Wigner limit" => "massless_wigner_limit.md",
    ],
)

deploydocs(;
    repo = "github.com/RUB-EP1/InstructionalDecayTrees.jl",
    root = DOCS,
    target = "build",
    versions = ["v#.#.#", "dev" => "dev"],
)
