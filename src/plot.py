import matplotlib.pyplot as plt
# import logging
plt.set_loglevel("info")


def get_best_methods(df_ns, metric_col, ascending):
    # Identify the best method for each group based on f max
    best_methods = []
    for group, df_g in df_ns.groupby(level='group'):
        best_methods.append(df_g[metric_col].idxmax())

    # Sort the best methods based on f
    df_best_methods = df_ns[df_ns['cov'] > 0].loc[best_methods].sort_values(by=['is_baseline', metric_col],
                                                                              ascending=[True, ascending])

    # Re-sort thresholds
    df_best_methods = df_best_methods.set_index('tau', append=True).groupby(
        level='label', group_keys=False, sort=False).apply(lambda x: x.sort_index(level='tau', ascending=False))

    return df_best_methods


def plot_curves(out_file, df, groups, ascending, metric_col, x_col, y_col, x_label, y_label, x_lim=None, y_lim=None, label_colors=None):
    """
    Plot precision recall curves
    """
    cmap = plt.get_cmap('tab10')
    fig, ax = plt.subplots(figsize=(10, 10))
    for method, df_m in df.groupby('label', sort=False):
        best = df_m.loc[df_m[metric_col].idxmax()] if ascending is False else df_m.loc[df_m[metric_col].idxmin()]
        # df_m.sort_values(by=[y_col, x_col], ascending=[False, True], inplace=True)
        # print(df_m)
        # color = cmap.colors[groups.get_loc(best.name[1]) % len(cmap.colors)]
        color = cmap.colors[groups.get_loc(best.label) % len(cmap.colors)]
        if label_colors is not None and method in label_colors:
            color = label_colors[method]
        # Curve
        ax.plot(df_m[x_col], df_m[y_col], '--' if best['is_baseline'] else '-',
                label="{} ({}={:.2f} C={:.2f})".format(method, metric_col.upper(), best[metric_col], best['cov']),
                color=color)
        # F-max circle
        plt.plot(best[x_col], best[y_col], 'o', color=color)

    ax.legend()
    # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 15})
    # ax.legend(bbox_to_anchor=(0, 1), loc='lower left', prop={'size': 15}, ncol=int(math.sqrt(len(labels)))-1)
    ax.set_ylabel(y_label) #, fontsize=18)
    ax.set_xlabel(x_label) #, fontsize=18)
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close('all')