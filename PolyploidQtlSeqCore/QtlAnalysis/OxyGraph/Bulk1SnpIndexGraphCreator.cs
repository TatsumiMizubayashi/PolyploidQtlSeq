using OxyPlot;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Bulk1 SnpIndexグラフクリエーター
    /// </summary>
    internal class Bulk1SnpIndexGraphCreator : GraphCreator
    {
        private const string TITLE = "Bulk1 SNP-index";
        private const string Y_AXIS_TITLE = "SNP-index";
        private const string FILE_BASENAME = "Bulk1 SNP-index";

        /// <summary>
        /// Bulk1 SnpIndexグラフクリエーターを作成する。
        /// </summary>
        /// <param name="config">設定</param>
        public Bulk1SnpIndexGraphCreator(GraphConfig config)
            : base(config, TITLE, Y_AXIS_TITLE, FILE_BASENAME)
        {
        }

        protected override void AddAnnotations(PlotModel plotModel)
        {
            // アノテーションなし
        }

        protected override void AddSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var snpIndesSeries = data.CreateBulk1SnpIndexScatterSeries(chrName);
            plotModel.Series.Add(snpIndesSeries);

            var zeroSnpIndexSeries = data.CreateBulk1ZeroSnpIndexScatterSeries(chrName);
            plotModel.Series.Add(zeroSnpIndexSeries);

            var averageSnpIndexSeries = data.CreateBulk1AverageSnpIndexLineSeries(chrName);
            foreach (var series in averageSnpIndexSeries)
            {
                plotModel.Series.Add(series);
            }
        }
    }
}
