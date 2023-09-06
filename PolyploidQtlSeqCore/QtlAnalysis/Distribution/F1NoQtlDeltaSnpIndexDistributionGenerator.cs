namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// F1 NoQtlΔSNP-index分布ジェネレータ
    /// </summary>
    internal class F1NoQtlDeltaSnpIndexDistributionGenerator
    {
        private readonly NoQtlDistributionSettings _settings;
        private readonly BulkF1SamplesGenerator _bulkF1SamplesGenerator;

        /// <summary>
        /// F1 NoQtlΔSNP-index分布ジェネレータを作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        public F1NoQtlDeltaSnpIndexDistributionGenerator(NoQtlDistributionSettings settings)
        {
            _settings = settings;

            var f1AltCountPatterns = F1GroupAltCountPattern.Generate(_settings.Ploidy, _settings.Parent2PlexNumber);
            _bulkF1SamplesGenerator = new BulkF1SamplesGenerator(f1AltCountPatterns, _settings.Ploidy);
        }


        /// <summary>
        /// F1 NoQtlΔSNP-index分布を生成する。
        /// </summary>
        /// <param name="bulk1">Bulk1</param>
        /// <param name="bulk2">Bulk2</param>
        /// <returns>F1 NoQtlΔSNP-index分布</returns>
        public F1NoQtlDeltaSnpIndexDistribution Generate(Bulk1 bulk1, Bulk2 bulk2)
        {
            var deltaSnpIndexList = new List<DeltaSnpIndex>();

            while (deltaSnpIndexList.Count <= _settings.ReplicatesNumber.Value)
            {
                var bulk1Samples = _bulkF1SamplesGenerator.Generate(_settings.Bulk1Number, bulk1);
                var bulk2Samples = _bulkF1SamplesGenerator.Generate(_settings.Bulk2Number, bulk2);

                var bulk1SnpIndex = bulk1Samples.SnpIndex;
                var bulk2SnpIndex = bulk2Samples.SnpIndex;
                var deltaSnpIndex = bulk2SnpIndex.CalcDeltaSnpIndex(bulk1SnpIndex);

                deltaSnpIndexList.Add(deltaSnpIndex);
            }

            return new F1NoQtlDeltaSnpIndexDistribution(deltaSnpIndexList);
        }
    }
}
