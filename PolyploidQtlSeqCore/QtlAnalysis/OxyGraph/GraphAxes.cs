namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ軸
    /// </summary>
    internal class GraphAxes
    {
        private readonly XAxisConfig _xConfig;
        private readonly YAxisConfig _bulkYConfig;
        private readonly YAxisConfig _deltaYConfig;
        private readonly YAxisConfig _scoreYConfig;
        private readonly YAxisConfig _countYConfig;

        /// <summary>
        /// グラフ軸インスタンスを作成する。
        /// </summary>
        /// <param name="x">X軸範囲</param>
        /// <param name="snpIndex">SnpIndex Y軸設定</param>
        /// <param name="delta">ΔSnpIndex Y軸設定</param>
        /// <param name="score">WindowScore Y軸設定</param>
        /// <param name="count">QTL数 Y軸設定</param>
        public GraphAxes(XAxisConfig x, YAxisConfig snpIndex, YAxisConfig delta, YAxisConfig score, YAxisConfig count)
        {
            _xConfig = x;
            _bulkYConfig = snpIndex;
            _deltaYConfig = delta;
            _scoreYConfig = score;
            _countYConfig = count;
        }

        /// <summary>
        /// Y軸補正を行う。
        /// </summary>
        /// <param name="setting">設定</param>
        /// <returns>補正済みGraphAxes</returns>
        public GraphAxes CorrectYAxis(GraphSettings setting)
        {
            var snpIndexGraphConfig = CreateSnpIndexGraphConfig(setting);
            var deltaSnpIndexGraphConfig = CreateDeltaSnpIndexGraphConfig(setting);
            var scoreGraphConfig = CreateWindowScoreGraphConfig(setting);
            var qtlCountGraphConfig = CreateQtlCountGraphConfig(setting);

            var creators = new GraphCreator[]
            {
                new BulkSnpIndexGraphCreator(snpIndexGraphConfig),
                new DeltaSnpIndexGraphCreator(deltaSnpIndexGraphConfig),
                new WindowScoreGraphCreator(scoreGraphConfig),
                new QtlVariantCountGraphCreator(qtlCountGraphConfig)
            };

            var yAxisWidths = creators.Select(x => x.GetYAxisWidth()).ToArray();
            var maxWidth = yAxisWidths.MaxBy(x => x.Width) ?? throw new InvalidOperationException("maxWidthがnull");
            var correctYConfigs = yAxisWidths.Select(x => x.Correct(maxWidth)).ToArray();

            return new GraphAxes(
                _xConfig,
                correctYConfigs[0],
                correctYConfigs[1],
                correctYConfigs[2],
                correctYConfigs[3]);
        }

        /// <summary>
        /// Bulk SnpIndexグラフ設定を作成する。
        /// </summary>
        /// <param name="setting">設定</param>
        /// <returns>Bulk SnpIndexグラフ設定</returns>
        public GraphConfig CreateSnpIndexGraphConfig(GraphSettings setting)
        {
            return new GraphConfig(_xConfig, _bulkYConfig, setting);
        }

        /// <summary>
        /// ΔSnpIndexグラフ設定を作成する。
        /// </summary>
        /// <param name="setting">設定</param>
        /// <returns>ΔSnpIndexグラフ設定</returns>
        public GraphConfig CreateDeltaSnpIndexGraphConfig(GraphSettings setting)
        {
            return new GraphConfig(_xConfig, _deltaYConfig, setting);
        }


        /// <summary>
        /// WindowScoreグラフ設定を作成する。
        /// </summary>
        /// <param name="setting">設定</param>
        /// <returns>WindowScoreグラフ設定</returns>
        public GraphConfig CreateWindowScoreGraphConfig(GraphSettings setting)
        {
            return new GraphConfig(_xConfig, _scoreYConfig, setting);
        }

        /// <summary>
        /// QTL数グラフ設定を作成する。
        /// </summary>
        /// <param name="setting">設定</param>
        /// <returns>QTL数グラフ設定</returns>
        public GraphConfig CreateQtlCountGraphConfig(GraphSettings setting)
        {
            return new GraphConfig(_xConfig, _countYConfig, setting);
        }
    }
}
