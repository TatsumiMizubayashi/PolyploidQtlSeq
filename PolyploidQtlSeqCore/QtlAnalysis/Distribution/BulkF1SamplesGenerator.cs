namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// Bulk F1サンプルジェネレータ
    /// </summary>
    internal class BulkF1SamplesGenerator
    {
        private readonly int[] _f1AltCountPatterns;
        private readonly Ploidy _ploidy;
        private readonly RandomSampling _randomSampling;

        /// <summary>
        /// Bulk F1サンプルジェネレータを作成する。
        /// </summary>
        /// <param name="f1AltCountPatterns">F1集団Alt数パターン配列</param>
        /// <param name="ploidy">ploidy</param>
        public BulkF1SamplesGenerator(int[] f1AltCountPatterns, Ploidy ploidy)
        {
            _f1AltCountPatterns = f1AltCountPatterns;
            _ploidy = ploidy;
            _randomSampling = new RandomSampling();
        }

        /// <summary>
        /// Bulk1 F1サンプルを生成する。
        /// </summary>
        /// <param name="bulk1Number">Bulk1個体数</param>
        /// <param name="bulk1">Bulk1情報</param>
        /// <returns>Bulk1 F1サンプル</returns>
        public Bulk1F1Samples Generate(Bulk1Number bulk1Number, Bulk1 bulk1)
        {
            var f1Samples = GenerateF1Samples(bulk1Number.Value);

            return new Bulk1F1Samples(f1Samples, bulk1);
        }

        /// <summary>
        /// Bulk2 F1サンプルを生成する。
        /// </summary>
        /// <param name="bulk2Number">Bulk2個体数</param>
        /// <param name="bulk2">Bulk2情報</param>
        /// <returns>Bulk2 F1サンプル</returns>
        public Bulk2F1Samples Generate(Bulk2Number bulk2Number, Bulk2 bulk2)
        {
            var f1Samples = GenerateF1Samples(bulk2Number.Value);

            return new Bulk2F1Samples(f1Samples, bulk2);
        }

        /// <summary>
        /// F1集団パターンから個体数分ランダムサンプリングしたF1サンプルを生成する。
        /// </summary>
        /// <param name="bulkNumber">個体数</param>
        /// <returns>F1サンプル</returns>
        private F1Samples GenerateF1Samples(int bulkNumber)
        {
            var bulkSampleAltCountArray = _randomSampling.Choice(_f1AltCountPatterns, bulkNumber);

            return new F1Samples(bulkSampleAltCountArray, _ploidy);
        }
    }
}
