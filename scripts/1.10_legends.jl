let
    p = plot(;title = "ider_colors")
    for (label, color) in ider_colors
        scatter!(p, [0.0], [0.0]; color, label)
    end
    mysavefig(p, "ider_color")
end

let
    p = plot(;title = "source_markers")
    for (src, m) in source_markers
        label = string(nameof(src))
        scatter!(p, [0.0], [0.0]; 
            m, label, color = :black
        )
    end
    mysavefig(p, "source_markers")
end

let
    p = plot(;title = "method_colors")
    for (method, m) in method_colors
        label = string(method)
        scatter!(p, [0.0], [0.0]; 
            m, label, color = :black
        )
    end
    mysavefig(p, "method_colors")
end

