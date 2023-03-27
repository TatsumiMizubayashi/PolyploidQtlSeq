using OxyPlot;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Bulk SnpIndexグラフクリエーター
    /// </summary>
    internal class BulkSnpIndexGraphCreator : GraphCreator
    {
        private const string TITLE = "Bulk SNP-index";
        private const string Y_AXIS_TITLE = "SNP-index";
        private const string FILE_BASENAME = "Bulk SNP-index";

        /// <summary>
        /// Bulk SnpIndexグラフクリエーターを作成する。
        /// </summary>
        /// <param name="config">設定</param>
        public BulkSnpIndexGraphCreator(GraphConfig config)
            : base(config, TITLE, Y_AXIS_TITLE, FILE_BASENAME)
        {
        }

        protected override void AddAnnotations(PlotModel plotModel)
        {
            // アノテーションなし
        }

        protected override void AddSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            AddBulk2ScatterSeries(plotModel, chrName, data);
            AddBulk1ScatterSeries(plotModel, chrName, data);

            AddBulk2LineSeries(plotModel, chrName, data);
            AddBulk1LineSeries(plotModel, chrName, data);
        }

        /// <summary>
        /// Bulk1のScatterSeriesを追加する。
        /// </summary>
        /// <param name="plotModel">PlotModel</param>
        /// <param name="chrName">染色体名</param>
        /// <param name="data">グラフデータ</param>
        private static void AddBulk1ScatterSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var snpIndesSeries = data.CreateBulk1SnpIndexScatterSeries(chrName);
            plotModel.Series.Add(snpIndesSeries);

            var zeroSnpIndexSeries = data.CreateBulk1ZeroSnpIndexScatterSeries(chrName);
            plotModel.Series.Add(zeroSnpIndexSeries);
        }

        /// <summary>
        /// Bulk1のLineSeriesを追加する。
        /// </summary>
        /// <param name="plotModel">PlotModel</param>
        /// <param name="chrName">染色体名</param>
        /// <param name="data">グラフデータ</param>
        private static void AddBulk1LineSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var averageSnpIndexSeries = data.CreateBulk1AverageSnpIndexLineSeries(chrName);
            foreach (var series in averageSnpIndexSeries)
            {
                plotModel.Series.Add(series);
            }
        }

        /// <summary>
        /// Bulk2のシリーズを追加する。
        /// </summary>
        /// <param name="plotModel">PlotModel</param>
        /// <param name="chrName">染色体名</param>
        /// <param name="data">グラフデータ</param>
        private static void AddBulk2ScatterSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var snpIndesSeries = data.CreateBulk2SnpIndexScatterSeries(chrName);
            plotModel.Series.Add(snpIndesSeries);

            var zeroSnpIndexSeries = data.CreateBulk2ZeroSnpIndexScatterSeries(chrName);
            plotModel.Series.Add(zeroSnpIndexSeries);

        }

        /// <summary>
        /// Bulk2のLineSeriesを追加する。
        /// </summary>
        /// <param name="plotModel">PlotModel</param>
        /// <param name="chrName">染色体名</param>
        /// <param name="data">グラフデータ</param>
        private static void AddBulk2LineSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var averageSnpIndexSeries = data.CreateBulk2AverageSnpIndexLineSeries(chrName);
            foreach (var series in averageSnpIndexSeries)
            {
                plotModel.Series.Add(series);
            }
        }
    }
}
