namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// AD
    /// </summary>
    internal class AD
    {
        private static readonly char[] _delimiter = [','];
        private const string NO_DATA = ".";

        private readonly ReadCounter _readCounter;

        /// <summary>
        /// ADインスタンスを作成する。
        /// </summary>
        /// <param name="ad">ADの値</param>
        public AD(string ad)
        {
            ArgumentException.ThrowIfNullOrEmpty(ad);

            Value = ad;
            IsNoData = ad.Contains(NO_DATA);
            ReadCounts = IsNoData ? [] : [.. ad.Split(_delimiter).Select(x => int.Parse(x))];
            Depth = ReadCounts.Sum();
            _readCounter = ReadCounter.Create(VariantAlleleType.RefAlt);    // 一番頻度が多いタイプを既定とする
        }

        /// <summary>
        /// ADインスタンスを作成する。
        /// </summary>
        /// <param name="ad">AD</param>
        /// <param name="counter">リードカウンター</param>
        public AD(AD ad, ReadCounter counter)
        {
            Value = ad.Value;
            IsNoData = ad.IsNoData;
            ReadCounts = ad.ReadCounts;
            Depth = ad.Depth;
            _readCounter = counter;
        }

        /// <summary>
        /// ADの値を取得する。
        /// </summary>
        public string Value { get; }

        /// <summary>
        /// NoDataかどうかを取得する。
        /// </summary>
        public bool IsNoData { get; }

        /// <summary>
        /// 各アレルのリード数を取得する。
        /// </summary>
        public int[] ReadCounts { get; }

        /// <summary>
        /// Depthを取得する。
        /// </summary>
        public int Depth { get; }

        /// <summary>
        /// Alleleカウントを取得する。
        /// </summary>
        /// <returns>(refCount, altCount)</returns>
        public (int refCount, int altCount) GetAlleleCount()
        {
            return _readCounter.Count(this);
        }
    }
}
