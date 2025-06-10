import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import upset_plot


# Import your separate modules (uncomment when ready)
# import upset_plot_module
# import box_plot_module

def main():
    st.set_page_config(
        page_title="CMMC Analysis Dashboard",
        page_icon="🧬",
        layout="wide"
    )

    st.title("🧬 CMMC Analysis Dashboard")
    st.markdown("---")

    # Sidebar configuration
    with st.sidebar:
        st.header("📊 Analysis Configuration")

        # CMMC Enrichment Task ID input
        st.subheader("CMMC Enrichment")
        cmmc_task_id = st.text_input(
            "Task ID",
            placeholder="Enter CMMC Enrichment Task ID",
            help="Input your CMMC enrichment task identifier"
        )

        st.markdown("---")

        # FBMN Task ID input
        st.subheader("FBMN Analysis")
        fbmn_task_id = st.text_input(
            "FBMN Task ID",
            placeholder="Enter FBMN Task ID",
            help="Input your Feature-Based Molecular Network task identifier"
        )

        st.markdown("---")

        # Metadata table upload
        st.subheader("📁 Data Upload")
        uploaded_file = st.file_uploader(
            "Upload Metadata Table",
            type=['csv', 'xlsx', 'tsv'],
            help="Upload your metadata table (CSV, Excel, or TSV format)"
        )

        # Display upload status
        if uploaded_file is not None:
            st.success(f"✅ File uploaded: {uploaded_file.name}")
            try:
                # Read the uploaded file
                if uploaded_file.name.endswith('.csv'):
                    metadata_df = pd.read_csv(uploaded_file)
                elif uploaded_file.name.endswith('.tsv'):
                    metadata_df = pd.read_csv(uploaded_file, sep='\t')
                else:  # Excel files
                    metadata_df = pd.read_excel(uploaded_file)

                st.info(f"📋 Rows: {len(metadata_df)} | Columns: {len(metadata_df.columns)}")

                # Show preview
                with st.expander("Preview Data"):
                    st.dataframe(metadata_df.head(), use_container_width=True)

            except Exception as e:
                st.error(f"❌ Error reading file: {str(e)}")
        else:
            st.info("📤 Please upload a metadata table")

        st.markdown("---")

        # Analysis button
        run_analysis = st.button(
            "🚀 Run Analysis",
            type="primary",
            use_container_width=True,
            disabled=not (cmmc_task_id and fbmn_task_id and uploaded_file)
        )

    # Main content area
    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("📈 UpSet Plot")
        upset_fig = upset_plot.generate_upset_plot(pd.read_csv("~/Downloads/upset_data.tsv", sep='\t')
                                                   )
        st.pyplot(upset_fig, use_container_width=True)

        if run_analysis and uploaded_file is not None:
            # Placeholder for upset plot generation
            upset_fig = upset_plot.generate_upset_plot(pd.read_csv("~/Downloads/upset_data.tsv")
            )
            st.pyplot(upset_fig, use_container_width=True)

            # Temporary placeholder visualization
            st.info("🔄 UpSet plot will be generated here using your separate module")

            # Create a simple placeholder plot
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=[1, 2, 3, 4],
                y=[10, 11, 12, 13],
                mode='markers+lines',
                name='Placeholder'
            ))
            fig.update_layout(
                title="UpSet Plot Placeholder",
                xaxis_title="Intersections",
                yaxis_title="Count",
                height=400
            )
            st.plotly_chart(fig, use_container_width=True)

        else:
            st.info("📋 Configure parameters in sidebar and run analysis to display UpSet plot")

    with col2:
        st.subheader("📊 Box Plots")

        if run_analysis and uploaded_file is not None:
            # Placeholder for box plot generation
            # box_fig = box_plot_module.generate_box_plots(
            #     cmmc_task_id=cmmc_task_id,
            #     fbmn_task_id=fbmn_task_id,
            #     metadata=metadata_df
            # )
            # st.plotly_chart(box_fig, use_container_width=True)

            # Temporary placeholder visualization
            st.info("🔄 Box plots will be generated here using your separate module")

            # Create a simple placeholder plot
            fig = go.Figure()
            fig.add_trace(go.Box(
                y=[1, 2, 3, 4, 5, 6, 7, 8, 9],
                name="Sample 1"
            ))
            fig.add_trace(go.Box(
                y=[2, 3, 4, 5, 6, 7, 8, 9, 10],
                name="Sample 2"
            ))
            fig.update_layout(
                title="Box Plot Placeholder",
                yaxis_title="Values",
                height=400
            )
            st.plotly_chart(fig, use_container_width=True)

        else:
            st.info("📋 Configure parameters in sidebar and run analysis to display box plots")

    # Status information
    st.markdown("---")

    # Display current configuration
    with st.expander("🔧 Current Configuration", expanded=False):
        col1, col2, col3 = st.columns(3)

        with col1:
            st.write("**CMMC Task ID:**")
            st.code(cmmc_task_id if cmmc_task_id else "Not set")

        with col2:
            st.write("**FBMN Task ID:**")
            st.code(fbmn_task_id if fbmn_task_id else "Not set")

        with col3:
            st.write("**Metadata File:**")
            st.code(uploaded_file.name if uploaded_file else "Not uploaded")


if __name__ == "__main__":
    main()