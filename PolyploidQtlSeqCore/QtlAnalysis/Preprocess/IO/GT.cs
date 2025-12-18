namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// GT
    /// </summary>
    internal class GT
    {
        private static readonly char[] _delimiter = ['/', '|'];

        private const string ALLELE_DELIMITER = "/";
        private const string NO_DATA = ".";


        /// <summary>
        /// GTを作成する。
        /// </summary>
        /// <param name="gt">GT</param>
        public GT(string gt)
        {
            ArgumentException.ThrowIfNullOrEmpty(gt);

            Value = gt;
            IsNoData = gt.Contains(NO_DATA);
            Indexes = IsNoData
                ? []
                : [.. gt.Split(_delimiter).Select(x => int.Parse(x))];
            Type = GetGtType();
        }

        /// <summary>
        /// GTの値を取得する。
        /// </summary>
        public string Value { get; }

        /// <summary>
        /// NoDataかどうかを取得する。
        /// </summary>
        public bool IsNoData { get; }

        /// <summary>
        /// アレルindex配列を取得する。
        /// 重複は除外されている。
        /// </summary>
        public int[] Indexes { get; }

        /// <summary>
        /// GTの種類を取得する。
        /// </summary>
        public GtType Type { get; }

        /// <summary>
        /// アレル塩基を取得する。
        /// </summary>
        /// <param name="allAlleles">全アレル塩基配列</param>
        /// <returns>アレル塩基</returns>
        public string ToAllele(string[] allAlleles)
        {
            if (allAlleles.Length == 0) throw new ArgumentException(null, nameof(allAlleles));

            if (IsNoData) return NO_DATA;

            var uniqueIndexes = Indexes.Distinct().ToArray();
            if (uniqueIndexes.Length == 1) return allAlleles[uniqueIndexes[0]];

            var alleles = Indexes.Select(x => allAlleles[x]);
            return string.Join(ALLELE_DELIMITER, alleles);
        }

        /// <summary>
        /// GTの種類を取得する。
        /// </summary>
        /// <returns>GTの種類</returns>
        private GtType GetGtType()
        {
            if (IsNoData) return GtType.NoData;

            var uniqueIndexes = Indexes.Distinct().ToArray();
            if (uniqueIndexes.Length != 1) return GtType.Hetero;

            return uniqueIndexes[0] == 0 ? GtType.RefHomo : GtType.AltHomo;
        }
    }
}
