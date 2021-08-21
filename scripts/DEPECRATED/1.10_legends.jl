let
    fontsize = 13
    commomm_params = (;dpi = 1000,
        thickness_scaling = 1.3, 
        xguidefontsize = fontsize, 
        yguidefontsize = fontsize
    )

    # let
    #     p = plot(;title = "ider_colors", commomm_params...)
    #     for (label, color) in ider_colors
    #         scatter!(p, [0.0], [0.0]; color, label)
    #     end
    #     mysavefig(p, "ider_color")
    # end

    # let
    #     p = plot(;title = "source_markers", commomm_params...)
    #     srcs = [Kd, Nd, Fd] # cronological order
    #     for src in srcs
    #         marker = source_markers[src]
    #         label = source_labels[src]
    #         scatter!(p, [0.0], [0.0]; 
    #             marker, label, color = :black
    #         )
    #     end
    #     mysavefig(p, "source_markers")
    # end

    # let
    #     p = plot(;title = "subs_markers", commomm_params...)
    #     for (sub, m) in subs_markers
    #         label = string(sub)
    #         scatter!(p, [0.0], [0.0]; 
    #             m, label, color = :black
    #         )
    #     end
    #     mysavefig(p, "subs_markers")
    # end

    # let
    #     p = plot(;title = "method_colors", commomm_params...)
    #     for (method, m) in method_colors
    #         label = string(method)
    #         scatter!(p, [0.0], [0.0]; 
    #             m, label, color = :black
    #         )
    #     end
    #     mysavefig(p, "method_colors")
    # end

    # let
    #     p = plot(;title = "method_colors", commomm_params...)
    #     for (method, m) in method_colors
    #         label = string(method)
    #         scatter!(p, [0.0], [0.0]; 
    #             m, label, color = :black
    #         )
    #     end
    #     mysavefig(p, "method_colors")
    # end

    let
        exp_label = "\$ \\textbf{exp.} \$"
        me_label = "\$\\textbf{ME}^{\\langle \\textbf{1} \\rangle}\$"
        fba_label = "\$ \\textbf{FBA} \$"
        p = plot(;title = "method_comparizon", commomm_params...)
        for (label, marker) in [
            (exp_label, (:star, :black)),
            (fba_label, (:square, :red)),
            (me_label, (:circle, :blue)),
        ]
        scatter!(p, [0.0], [0.0]; label, marker, 
            grid = false, 
        )
        vspan!(p, [-1.0, 1.0]; label = "",
            linecolor = :darkgray, fillcolor = :darkgray, 
            fillalpha = 0.6, linealpha = 0.6, 
            xlim = [-1, 1],
            ylim = [-1, 1],
        )
        end

        mysavefig(p, "method_comparizon")
    end
end