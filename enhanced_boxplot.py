from itertools import combinations
import io
import zipfile

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy.stats import mannwhitneyu, kruskal, ttest_ind, f_oneway
from statsmodels.stats.multitest import multipletests

from utils import insert_contribute_link, insert_request_dep_correction_link, render_details_card

ORIGIN_LIST = [
    "Ambiguous",
    "De novo biosynthesis by microbes",
    "Diet",
    "Drug",
    "Exposure",
    "Exposure/diet",
    "Host",
    "Host metabolism of microbial metabolites",
    "Insecticides/pesticides",
    "Microbial metabolism of drugs",
    "Microbial metabolism of food molecules",
    "Microbial metabolism of host-derived molecules",
    "Microbial metabolism of microbial-derived molecules",
    "Microbial metabolism of other human-made molecules",
    "Unknown/Undefined",
]
SOURCE_LIST = [
    "Microbial",
    "Host",
    "Diet",
    "Unknown",
    "Ambiguous",
    "Drug",
    "Exposure",
    "Pesticides/insecticides",
    "Other human-made molecules",
]


def render_plot_style_options(groups, color_prefix="color", rotation_key="labels_rot"):
    check_col, logscale_col = st.columns(2)
    with check_col:
        custom_check = st.checkbox(
            "Use custom colors", key='stats_custom_colors',
        )
    with logscale_col:
        logscale_check = st.checkbox(
            "Use log scale for y-axis",
            key='stats_log_scale',
            help="Enable to use logarithmic scale for y-axis (useful for skewed distributions)",
        )
    if custom_check:
        color_cols = st.columns(3)
        for idx, item in enumerate(groups):
            with color_cols[idx % 3]:
                st.color_picker(
                    f"{item}",
                    key=f"{color_prefix}_{item}",
                    help="Select a color for the group",
                    value=st.session_state.get(f"{color_prefix}_{item}", "#1f77b4"),
                )
    rotate_labels_angle = st.slider(
        "Rotate x-axis labels",
        min_value=-90,
        max_value=90,
        value=0,
        step=45,
        key=rotation_key,
    )
    return custom_check, logscale_check, rotate_labels_angle


def perform_statistical_test(feature_data, grouping_column, intensity_col, selected_groups, test_type, alpha=0.05,
                             stratify_column=None, selected_strata=None):
    """
    Perform statistical tests on grouped data with optional stratification
    """
    results = {}

    if stratify_column and selected_strata:
        # Perform tests for each stratum
        strata_results = {}

        for stratum in selected_strata:
            stratum_data = feature_data[feature_data[stratify_column] == stratum]

            # Get data for each group within this stratum
            group_data = []
            group_names = []

            for group in selected_groups:
                group_subset = stratum_data[stratum_data[grouping_column] == group][intensity_col].dropna()
                if len(group_subset) > 0:
                    group_data.append(group_subset)
                    group_names.append(group)

            if len(group_data) < 2:
                strata_results[stratum] = {"error": f"Need at least 2 groups with data for stratum: {stratum}"}
                continue

            try:
                stratum_result = _perform_single_test(group_data, group_names, test_type, alpha)
                stratum_result["stratum"] = stratum
                stratum_result["n_samples"] = {name: len(data) for name, data in zip(group_names, group_data)}
                strata_results[stratum] = stratum_result
            except Exception as e:
                strata_results[stratum] = {"error": f"Statistical test failed for {stratum}: {str(e)}"}

        results["stratified_results"] = strata_results
    else:
        # Original non-stratified analysis
        group_data = []
        group_names = []

        for group in selected_groups:
            group_subset = feature_data[feature_data[grouping_column] == group][intensity_col].dropna()
            if len(group_subset) > 0:
                group_data.append(group_subset)
                group_names.append(group)

        if len(group_data) < 2:
            return {"error": "Need at least 2 groups with data for statistical testing"}

        try:
            results = _perform_single_test(group_data, group_names, test_type, alpha)
        except Exception as e:
            return {"error": f"Statistical test failed: {str(e)}"}

    return results


def _perform_single_test(group_data, group_names, test_type, alpha):
    """Helper function to perform a single statistical test"""
    results = {}

    if test_type == "Mann-Whitney U (2 groups)":
        if len(group_data) != 2:
            return {"error": "Mann-Whitney U test requires exactly 2 groups"}
        statistic, p_value = mannwhitneyu(group_data[0], group_data[1], alternative='two-sided')
        results = {
            "test": "Mann-Whitney U",
            "statistic": statistic,
            "p_value": p_value,
            "significant": p_value < alpha,
            "groups": group_names[:2]
        }

    elif test_type == "Kruskal-Wallis (>= 3 groups)":
        statistic, p_value = kruskal(*group_data)
        results = {
            "test": "Kruskal-Wallis",
            "statistic": statistic,
            "p_value": p_value,
            "significant": p_value < alpha,
            "groups": group_names
        }

    elif test_type == "T-test (2 groups)":
        if len(group_data) != 2:
            return {"error": "T-test requires exactly 2 groups"}
        statistic, p_value = ttest_ind(group_data[0], group_data[1])
        results = {
            "test": "Independent T-test",
            "statistic": statistic,
            "p_value": p_value,
            "significant": p_value < alpha,
            "groups": group_names[:2]
        }

    elif test_type == "ANOVA (>= 3 groups)":
        statistic, p_value = f_oneway(*group_data)
        results = {
            "test": "One-way ANOVA",
            "statistic": statistic,
            "p_value": p_value,
            "significant": p_value < alpha,
            "groups": group_names
        }

    # Add pairwise comparisons for multiple group tests
    if len(group_data) > 2 and test_type in ["Kruskal-Wallis (>= 3 groups)", "ANOVA (>= 3 groups)"]:
        pairwise_results = []
        for i, j in combinations(range(len(group_data)), 2):
            if test_type == "Kruskal-Wallis (>= 3 groups)":
                stat, p_val = mannwhitneyu(group_data[i], group_data[j], alternative='two-sided')
            else:
                stat, p_val = ttest_ind(group_data[i], group_data[j])

            pairwise_results.append({
                "group1": group_names[i],
                "group2": group_names[j],
                "p_value": p_val,
                "significant": p_val < alpha
            })

        # Apply multiple testing correction
        p_values = [result["p_value"] for result in pairwise_results]
        rejected, p_corrected, _, _ = multipletests(p_values, alpha=alpha, method='bonferroni')

        for i, result in enumerate(pairwise_results):
            result["p_corrected"] = p_corrected[i]
            result["significant_corrected"] = rejected[i]

        results["pairwise_comparisons"] = pairwise_results

    return results


def create_stratified_boxplot(df, feature_id, grouping_column, selected_groups, stratify_column=None,
                              selected_strata=None, test_results=None, color_mapping=None):
    """
    Create stratified boxplots with paired side-by-side layout in a single figure
    """
    # Filter data for selected groups
    plot_data = df[df[grouping_column].isin(selected_groups)].copy()

    if stratify_column and selected_strata:
        plot_data = plot_data[plot_data[stratify_column].isin(selected_strata)]

    # Get intensity data for the feature
    intensity_col = "Peak Area"  # Adjust based on your data structure
    if intensity_col not in plot_data.columns:
        possible_cols = [col for col in plot_data.columns if 'area' in col.lower() or 'intensity' in col.lower()]
        if possible_cols:
            intensity_col = possible_cols[0]
        else:
            st.error("Could not find intensity/peak area column")
            return None, None

    # Calculate global y-axis range for consistency
    y_min = plot_data[intensity_col].min()
    y_max = plot_data[intensity_col].max()
    y_range = y_max - y_min
    y_padding = y_range * 0.05  # 5% padding

    # Create single figure with all boxplots
    fig = go.Figure()

    # Default colors
    colors = color_mapping if color_mapping else px.colors.qualitative.Set1

    if stratify_column and selected_strata:
        fig = px.box(
            plot_data,
            x=grouping_column,
            y=intensity_col,
            category_orders={grouping_column: selected_groups,
                             stratify_column: selected_strata},
            color=grouping_column,
            color_discrete_sequence=px.colors.qualitative.Set1,
            color_discrete_map=None if not color_mapping else color_mapping,
            facet_col=stratify_column,
            points="all",
            hover_name='filename',
        )
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        facet_values = plot_data[stratify_column].dropna().unique()
        for i, value in enumerate(facet_values, 1):
            fig.update_yaxes(showticklabels=True, col=i)

        fig.update_layout(
            showlegend=False,
            xaxis_title=grouping_column,
            yaxis_title="Peak Area",
            boxgroupgap=0,
        )

    else:
        # Original single plot (no stratification)
        for i, group in enumerate(selected_groups):
            group_data = plot_data[plot_data[grouping_column] == group]

            color = colors.get(group) if isinstance(colors, dict) else colors[i % len(colors)]

            fig.add_trace(go.Box(
                y=group_data[intensity_col],
                name=group,
                boxpoints='all',
                pointpos=0,
                marker_color=color,
                line=dict(width=2),
                opacity=0.8
            ))

        # Update layout for single plot
        fig.update_layout(
            title=f"Statistical Comparison: {plot_data[plot_data['featureID'] == feature_id]['input_name'].iloc[0] if 'input_name' in plot_data.columns else f'Feature {feature_id}'}",
            xaxis_title=grouping_column,
            yaxis_title="Peak Area" if intensity_col else "Intensity",

            showlegend=False,
            height=600,
            template="plotly_white"
        )

        # Add statistical annotations for single plot
        if test_results and "error" not in test_results:
            p_value = test_results.get("p_value", 0)
            significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"

            fig.add_annotation(
                x=0.5, y=y_max + y_range * 0.1,
                xref="paper", yref="y",
                text=f"{test_results.get('test', 'Test')}: p = {p_value:.4f} ({significance})",
                showarrow=False,
                font=dict(size=12),
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="black",
                borderwidth=1
            )

    fig.update_traces(
        boxpoints='all',
        marker=dict(opacity=0.6),
        pointpos=0,
        jitter=0.3,
    )
    fig.update_layout(yaxis=dict(exponentformat="power", showexponent="all"))

    return fig, plot_data


def add_pair_annotations(fig, plot_data, selected_strata, intensity_col, stratify_column, stats_by_facet):
    # facet columns (left→right)
    facet_vals = list(selected_strata)
    stats_by_facet = {k: [v] for k, v in stats_by_facet.items()}
    for col_idx, facet_val in enumerate(facet_vals, start=1):
        df_f = plot_data[plot_data[stratify_column] == facet_val]
        if df_f.empty or facet_val not in stats_by_facet:
            continue

        y_max = float(df_f[intensity_col].max())
        step = 0.05 * (df_f[intensity_col].max() - df_f[intensity_col].min() + 1e-9)
        y_level = y_max + step
        if not stats_by_facet[facet_val][0].get("error", False):
            for i, d in enumerate(stats_by_facet[facet_val]):
                g1, g2 = d["groups"][0], d["groups"][1]
                text = f'{d["test"]}: p={d["p_value"]:.3g}'
                y = y_level + i * step  # stack multiple pairs
                # annotation text (place near the right group)
                fig.add_annotation(
                    x=g2, y=y + step * 0.2, text=text,
                    showarrow=False, xanchor="center", yanchor="bottom",
                    row=1, col=col_idx
                )

    # extend y-limits a bit so annotations fit
    fig.update_yaxes(matches=None)  # allow independent range per facet
    for col_idx in range(1, len(facet_vals) + 1):
        fig.update_yaxes(row=1, col=col_idx, autorange=True)


def generate_all_feature_plots_zip(filtered_df, fid_items, grouping_column, selected_groups,
                                 stratify_column, selected_strata, selected_test, alpha_level,
                                 custom_colors=None, use_log_scale=False, rotate_angle=0):
    """
    Generate plots for all features and return as a ZIP file in memory
    """
    try:
        zip_buffer = io.BytesIO()

        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            progress_bar = st.progress(0)

            for idx, feature_item in enumerate(fid_items):
                try:
                    # Update progress
                    progress = (idx + 1) / len(fid_items)
                    progress_bar.progress(progress, text=f"Generating plot {idx + 1}/{len(fid_items)}: {feature_item[:50]}...")

                    # Extract feature ID
                    feature_id = int(feature_item.split(":")[0]) if ":" in feature_item else feature_item

                    # Filter data for this feature
                    filter_conditions = [
                        (filtered_df['featureID'] == feature_id),
                        (filtered_df[grouping_column].isin(selected_groups))
                    ]
                    if stratify_column and selected_strata:
                        filter_conditions.append(filtered_df[stratify_column].isin(selected_strata))

                    feature_data = filtered_df[
                        filter_conditions[0] & filter_conditions[1] &
                        (filter_conditions[2] if len(filter_conditions) > 2 else True)
                    ].copy()

                    if feature_data.empty:
                        continue  # Skip features with no data

                    # Find intensity column
                    intensity_col = 'Peak area'
                    if intensity_col not in feature_data.columns:
                        possible_cols = [col for col in feature_data.columns if 'area' in col.lower() or 'intensity' in col.lower()]
                        if possible_cols:
                            intensity_col = possible_cols[0]
                        else:
                            continue  # Skip if no intensity column found

                    # Apply log scale transformation if needed
                    if use_log_scale:
                        feature_data[intensity_col] = feature_data[intensity_col].apply(
                            lambda x: x if x > 0 else 1e-9  # Avoid log(0)
                        )

                    # Perform statistical test
                    test_results = perform_statistical_test(
                        feature_data, grouping_column, intensity_col, selected_groups,
                        selected_test, alpha_level, stratify_column, selected_strata
                    )

                    # Create plot
                    fig, plot_data = create_stratified_boxplot(
                        feature_data, feature_id, grouping_column, selected_groups,
                        stratify_column, selected_strata, test_results, custom_colors
                    )

                    if fig:
                        # Apply styling
                        fig.update_xaxes(tickangle=rotate_angle)
                        if use_log_scale:
                            fig.update_yaxes(type="log", exponentformat="power", showexponent="all")

                        # Add annotations for stratified plots
                        if stratify_column:
                            add_pair_annotations(fig, plot_data, selected_strata, intensity_col,
                                               stratify_column, test_results.get("stratified_results", {}))

                        # Convert to SVG
                        svg_bytes = fig.to_image(format="svg")

                        # Create filename
                        safe_feature_name = str(feature_item).replace(":", "_").replace("/", "_").replace("\\", "_")[:50]
                        filename = f"plot_{idx+1:03d}_{feature_id}_{safe_feature_name}.svg"

                        # Add to ZIP
                        zip_file.writestr(filename, svg_bytes)

                except Exception as e:
                    # Log error but continue with other features
                    st.warning(f"Error generating plot for feature {feature_item}: {str(e)}")
                    continue

            progress_bar.empty()

        zip_buffer.seek(0)
        return zip_buffer.getvalue()

    except Exception as e:
        st.error(f"Error creating ZIP file: {str(e)}")
        return None


def render_statistical_boxplot_tab(merged_df):
    """
    Render the enhanced statistical boxplot tab with stratification
    """
    st.subheader("📊 Statistical Boxplot Analysis",
                 help="Perform statistical comparisons between groups with optional stratification and customizable tests")

    if merged_df is None or merged_df.empty:
        st.warning("No data available. Please run the analysis first.")
        return

    # Configuration columns
    config_col1, config_col2, config_col3, config_col4 = st.columns(4)

    with config_col1:
        # Select grouping variable
        metadata_columns = [col for col in merged_df.columns if
                            col.startswith("ATTRIBUTE_")]
        grouping_column = st.selectbox(
            ":blue-badge[Step 1] Primary Grouping",
            metadata_columns,
            help="Select the primary column to group samples by"
        )

    with config_col2:
        # Select groups to compare
        if grouping_column:
            available_groups = sorted(merged_df[grouping_column].dropna().unique())
            selected_groups = st.multiselect(
                ":orange-badge[Step 2] Groups to Compare",
                available_groups,
                default=available_groups[:min(2, len(available_groups))],
                help="Select 2 or more groups for statistical comparison"
            )

    with config_col3:
        # Select stratification variable
        stratify_options = [None] + [col for col in metadata_columns if col != grouping_column]
        stratify_column = st.selectbox(
            ":violet-badge[Step 4A] Stratify by (Optional)",
            stratify_options,
            format_func=lambda x: "None" if x is None else x,
            help="Optional: Select a column to create paired boxplots for each category in a single figure"
        )

    with config_col4:
        # Select strata if stratification is enabled
        selected_strata = []
        if stratify_column:
            available_strata = sorted(merged_df[stratify_column].dropna().unique())
            selected_strata = st.multiselect(
                ":violet-badge[Step 4B] Categories to Include",
                available_strata,
                default=available_strata[:min(4, len(available_strata))],
                help="Select categories for paired side-by-side analysis"
            )

    # Statistical test selection
    test_col, alpha_col, filter_col = st.columns([1, 1, 2])

    with test_col:
        test_options = [
            "Mann-Whitney U (2 groups)",
            "Kruskal-Wallis (>= 3 groups)",
            "T-test (2 groups)",
            "ANOVA (>= 3 groups)"
        ]

        selected_test = st.selectbox(
            ":red-badge[Step 3A] Statistical Test",
            test_options,
            help="Choose appropriate test based on your data distribution and number of groups"
        )

    with alpha_col:
        alpha_level = st.slider(
            ":red-badge[Step 3B] Significance Level (α)",
            min_value=0.01,
            max_value=0.10,
            value=0.05,
            step=0.01,
            help="P-value threshold for statistical significance"
        )

    with filter_col:
        st.write(
            '<div style="height: 25px;"></div>', unsafe_allow_html=True
        )
        with st.expander(":gray-badge[Step 5] Filter by Source or Origin :gray-background[(Optional)]", icon="🔍"):
            # Add filtering options
            if 'input_molecule_origin' in merged_df.columns:
                origin_filter = st.multiselect(
                    "Molecule Origin",
                    ORIGIN_LIST,
                    help="Filter by molecule origin"
                )
            else:
                origin_filter = []

            if 'input_source' in merged_df.columns:
                source_filter = st.multiselect(
                    "Molecule Source",
                    SOURCE_LIST,
                    help="Filter by molecule source"
                )
            else:
                source_filter = []

    # Style and filter options
    style_col, feature_col = st.columns(2)

    with style_col:
        st.write(
            '<div style="height: 25px;"></div>', unsafe_allow_html=True
        )
        with st.expander("Style Options", icon=":material/palette:"):
            color_picker_prefix = "stats"
            use_custom_colors, use_log_scale, rotate_angle = render_plot_style_options(selected_groups, color_picker_prefix)
            custom_colors = {
                group: st.session_state.get(f"{color_picker_prefix}_{group}", "#1f77b4")
                for i, group in enumerate(selected_groups)
            }

    # Apply filters to data
    filtered_df = merged_df.copy()
    filter_string = " "
    if origin_filter:
        filtered_df = filtered_df[filtered_df['input_molecule_origin'].isin(origin_filter)]
        filter_string += f"Origin: {', '.join(origin_filter)}; "
    if source_filter:
        filtered_df = filtered_df[filtered_df['input_source'].isin(source_filter)]
        filter_string += f"Source: {', '.join(source_filter)}"

    with feature_col:
        # Create feature dictionary
        if 'input_name' in filtered_df.columns:
            feat_id_dict = dict(zip(filtered_df["featureID"], filtered_df["input_name"]))
            feat_id_dict = {str(k): str(v) for k, v in feat_id_dict.items()}
            feat_id_dict = dict(sorted(feat_id_dict.items(), key=lambda item: item[1]))
            fid_items = [f"{k}: {v}" for k, v in feat_id_dict.items()]
        else:
            feat_id_dict = {str(fid): f"Feature {fid}" for fid in filtered_df["featureID"].unique()}
            fid_items = list(feat_id_dict.values())

        selected_feature = st.selectbox(
            f":green-badge[Step 6] Feature to Plot :blue-badge[{len(fid_items)} features]" + f":red-badge[{filter_string}]",
            fid_items,
            index=0,
            help="Select the feature/metabolite to perform statistical analysis on"
        )

        # Add bulk download button
        if len(fid_items) > 1:
            if st.button(
                f":material/download: Download all plots (.zip)",
                help="Generate and download SVG plots for all features in the current selection",
                type="secondary",
                use_container_width=True
            ):
                with st.spinner(f"Generating {len(fid_items)} plots..."):
                    zip_data = generate_all_feature_plots_zip(
                        filtered_df, fid_items, grouping_column, selected_groups,
                        stratify_column, selected_strata, selected_test, alpha_level,
                        custom_colors if use_custom_colors else None, use_log_scale, rotate_angle
                    )

                if zip_data:
                    st.download_button(
                        label=":material/file_download: Download ZIP File",
                        data=zip_data,
                        file_name=f"all_feature_plots_{grouping_column}_{'paired' if stratify_column else 'simple'}.zip",
                        mime="application/zip",
                        type="primary"
                    )
                    st.success(f"Generated ZIP file with {len(fid_items)} plots successfully")

    # Main analysis section
    if selected_feature and selected_groups and len(selected_groups) >= 2:
        # Check stratification requirements
        if stratify_column and not selected_strata:
            st.warning("Please select categories to include for stratified analysis.")
            return

        feature_id = int(selected_feature.split(":")[0]) if ":" in selected_feature else selected_feature

        # Filter data for selected feature and groups
        filter_conditions = [
            (filtered_df['featureID'] == feature_id),
            (filtered_df[grouping_column].isin(selected_groups))
        ]

        if stratify_column and selected_strata:
            filter_conditions.append(filtered_df[stratify_column].isin(selected_strata))

        feature_data = filtered_df[
            filter_conditions[0] & filter_conditions[1] &
            (filter_conditions[2] if len(filter_conditions) > 2 else True)
            ].copy()

        if feature_data.empty:
            st.warning("No data found for the selected feature and groups combination.")
            return

        # Find the intensity column
        intensity_col = 'Peak area'
        if intensity_col not in feature_data.columns:
            possible_cols = [col for col in feature_data.columns if 'area' in col.lower() or 'intensity' in col.lower()]
            if possible_cols:
                intensity_col = possible_cols[0]
            else:
                st.error("Could not find intensity/peak area column")
                return

        # Perform statistical test
        with st.spinner("Performing statistical analysis..."):
            test_results = perform_statistical_test(
                feature_data, grouping_column, intensity_col, selected_groups,
                selected_test, alpha_level, stratify_column, selected_strata
            )

        plot_col, details_col = st.columns([3, 1])
        with plot_col:
            # sum to all peak areas if log scale is selected
            if use_log_scale:
                feature_data[intensity_col] = feature_data[intensity_col].apply(
                    lambda x: x if x > 0 else 1e-9  # Avoid log(0)
                )
            # Create and display the plot
            color_mapping = custom_colors if use_custom_colors else None
            fig, plot_data = create_stratified_boxplot(
                feature_data, feature_id, grouping_column, selected_groups,
                stratify_column, selected_strata, test_results, color_mapping
            )
            fig.update_xaxes(tickangle=rotate_angle)
            #apply log scale if selected
            if use_log_scale:
                fig.update_yaxes(type="log", exponentformat="power", showexponent="all")

            if fig:
                if stratify_column:
                    add_pair_annotations(fig, plot_data, selected_strata, intensity_col, stratify_column,
                                         test_results.get("stratified_results", {}))
                st.plotly_chart(fig, use_container_width=True)

                # Download options
                svg_bytes = fig.to_image(format="svg")
                col1, col2 = st.columns(2)
                with col1:
                    st.download_button(
                        "Download Plot (SVG)",
                        data=svg_bytes,
                        file_name=f"statistical_boxplot_{feature_id}_{grouping_column}_{'paired' if stratify_column else 'simple'}.svg",
                        mime="image/svg+xml",
                        icon=":material/download:"
                    )
                with col2:
                    csv_data = plot_data.to_csv(index=False)
                    st.download_button(
                        "Download Data (CSV)",
                        data=csv_data,
                        file_name=f"statistical_data_{feature_id}_{grouping_column}_{'paired' if stratify_column else 'simple'}.csv",
                        mime="text/csv",
                        icon=":material/download:"
                    )

        # Display statistical results in expander below the plot
        with st.expander(f"📊 Statistical Test Results - :blue-badge[{selected_test}]", expanded=True):
                if "stratified_results" in test_results:
                    # Display stratified results in columns
                    strata_items = list(test_results["stratified_results"].items())
                    n_cols = min(4, len(strata_items))
                    cols = st.columns(n_cols)

                    for i, (stratum, stratum_results) in enumerate(strata_items):
                        with cols[i % n_cols]:
                            with st.container(border=True):
                                st.write(f"**{stratify_column}: {stratum}**")

                                if "error" in stratum_results:
                                    st.error(stratum_results["error"])
                                else:
                                    st.metric(
                                        label="p-value",
                                        value=f"{stratum_results.get('p_value', 0):.4f}",
                                        delta="Sig" if stratum_results.get('significant', False) else "NS"
                                    )
                                    st.metric(
                                        label="Statistic",
                                        value=f"{stratum_results.get('statistic', 0):.3f}"
                                    )

                                    # Sample sizes for this stratum
                                    if "n_samples" in stratum_results:
                                        st.caption("Sample sizes:")
                                        for group, n in stratum_results["n_samples"].items():
                                            st.text(f"  {group}: n={n}")

                else:
                    # Display regular results in columns
                    col1, col2, col3 = st.columns(3)

                    if "error" in test_results:
                        with col1:
                            st.error(test_results["error"])
                    else:
                        with col1:
                            st.metric(
                                label=f"{test_results.get('test', 'Test')} p-value",
                                value=f"{test_results.get('p_value', 0):.4f}",
                                delta="Significant" if test_results.get('significant', False) else "Not significant"
                            )

                        with col2:
                            st.metric(
                                label="Test Statistic",
                                value=f"{test_results.get('statistic', 0):.4f}"
                            )

                        # Feature details in third column
                        with col3:
                            if 'input_name' in feature_data.columns:
                                st.subheader("Feature Details")
                                feature_info = feature_data.iloc[0]
                                for col in ['input_name', 'input_molecule_origin', 'input_source']:
                                    if col in feature_info:
                                        st.text(
                                            f"{col.replace('input_', '').replace('_', ' ').title()}: {feature_info[col]}")

        with details_col:
            # Show details card for the selected feature ID
            enriched_result = st.session_state.get("enriched_result")

            with st.expander("Details", icon=":material/info:"):
                columns_to_show = st.multiselect(
                    "Select columns to show in details card",
                    enriched_result.columns.tolist(),
                    default=[
                        "input_name",
                        "input_molecule_origin",
                        "input_source",
                    ],
                )
            render_details_card(
                enriched_result, int(feature_id), columns_to_show
            )
        insert_contribute_link(enriched_result, selected_feature)

        # renders a link to request a correction
        insert_request_dep_correction_link(enriched_result, selected_feature)

    else:
        st.info("Please select a feature and at least 2 groups to perform statistical analysis.")

    # Help section
    with st.expander("ℹ️ Help & Guide"):
        st.markdown("""
        **Enhanced Paired Boxplot Feature:**
        - Select a "Stratify by" column to create paired boxplots side by side in a single figure
        - All groups are shown for each stratum category, making comparisons easier
        - Useful for exploring interactions between different variables in one view
        - Statistical annotations appear above each stratum group
        
        **Statistical Test Selection Guide:**
        - **Mann-Whitney U (2 groups)**: Non-parametric test for comparing 2 independent groups
        - **Kruskal-Wallis (>= 3 groups)**: Non-parametric test for comparing 3+ independent groups
        - **T-test (2 groups)**: Parametric test for comparing 2 independent groups (assumes normality)
        - **ANOVA (>= 3 groups)**: Parametric test for comparing 3+ independent groups (assumes normality)
        
        **Significance Levels:**
        - *** p < 0.001, ** p < 0.01, * p < 0.05, ns: not significant

        """)


if __name__ == '__main__':
    file = 'merged_df_test.csv'
    test_options = [
        "Mann-Whitney U (2 groups)",
        "Kruskal-Wallis (>= 3 groups)",
        "T-test (2 groups)",
        "ANOVA (>= 3 groups)"
    ]

    merged_df = pd.read_csv(file)
    filtered_df = merged_df.copy()  # Assuming no filters for this test

    feature_id = 560
    # feature_id = st.selectbox('Feature ID', filtered_df['featureID'].unique(), help="Select a feature ID to analyze")

    feature_data = filtered_df[filtered_df['featureID'] == feature_id]  # Example feature ID

    # grouping_column = 'ATTRIBUTE_GF_SPF'
    grouping_column = st.selectbox("Grouping column",
                                   [col for col in feature_data.columns if "attribute" in col.lower()])

    # selected_groups = ['SPF', 'GF']
    available_groups = feature_data[grouping_column].unique()
    selected_groups = st.multiselect('Groups', available_groups,
                                     default=['SPF', 'GF'] if "SPF" in grouping_column else available_groups, )

    # intensity_col = 'Peak Area'
    intensity_col = st.selectbox("Intensity column", [col for col in feature_data.columns if
                                                      'area' in col.lower() or 'intensity' in col.lower()],
                                 help="Select the intensity/peak area column")

    selected_test = test_options[1]
    # selected_test = st.selectbox("Statistical Test", test_options, help="Choose the statistical test")

    # alpha_level = 0.05
    alpha_level = st.slider("Significance Level (α)", min_value=0.01, max_value=0.10, value=0.05, step=0.01,
                            help="P-value threshold for statistical significance")

    stratify_column = 'ATTRIBUTE_UBERONBodyPartName'
    # stratify_column = st.selectbox("Stratification column", [None] + [col for col in feature_data.columns if
    #                                                                   "attribute" in col.lower() and col != grouping_column],
    #                                format_func=lambda x: "None" if x is None else x,
    #                                index=0)

    # selected_strata = ['Ileum', 'Cecum', "Vagina"]
    selected_strata = ['Adrenal_gland', 'Bedding', "Bladder"]
    # if stratify_column:
    #     available_strata = feature_data[stratify_column].unique()
    #     selected_strata = st.multiselect("Body Part", available_strata, default=['Ileum', 'Cecum',
    #                                                                              "Stool"] if "BodyPart" in stratify_column else available_strata, )
    # else:
    #     selected_strata = []

    # use_custom_colors = False
    use_custom_colors = st.checkbox("Use custom colors", value=False)
    custom_colors = {}
    if use_custom_colors and selected_groups:
        for group in selected_groups:
            custom_colors[group] = st.color_picker(f"Color for {group}", value="#1f77b4")

    color_mapping = custom_colors if use_custom_colors else None

    test_results = perform_statistical_test(
        feature_data, grouping_column, intensity_col, selected_groups,
        selected_test, alpha_level, stratify_column, selected_strata
    )

    fig, plot_data = create_stratified_boxplot(
        feature_data, feature_id, grouping_column, selected_groups,
        stratify_column, selected_strata, test_results, color_mapping,
    )
    st.write(test_results)
    if selected_strata:
        add_pair_annotations(fig, plot_data, intensity_col, stratify_column, test_results.get("stratified_results", {}))

    st.plotly_chart(fig, use_container_width=True)

    st.write(merged_df.dtypes)
