"""Codex's plot grammar as helpers: 3x2 constrained grid, palette, faint grid, frameless legends, dpi 170."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PALETTE = dict(brown='#976b3b', teal='#147d72', blue='#546da8', grey='#252525', accent='#b03a2e')
BR, TE, BL, GY, AC = PALETTE['brown'], PALETTE['teal'], PALETTE['blue'], PALETTE['grey'], PALETTE['accent']
BY_N = {1: TE, 2: BL, 3: BR}


def grid(rows=3, cols=2, figsize=(12.5, 12)):
    return plt.subplots(rows, cols, figsize=figsize, layout='constrained')


def finish(fig, axes, title, legend_size=7):
    for a in axes.flat:
        a.grid(alpha=.2)
        if a.get_legend_handles_labels()[0]: a.legend(frameon=False, fontsize=legend_size)
    fig.suptitle(title)
    return fig
