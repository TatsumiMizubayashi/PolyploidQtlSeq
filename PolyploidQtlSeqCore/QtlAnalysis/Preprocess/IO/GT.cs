namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// GT
    /// </summary>
    internal class GT
    {
        private static readonly char[] _delimiter = ['/', '|'];

        private const string ALLELE_DELIMITER = "/";
        private const string NO_DATA_GT = "./.";
        private const string NO_DATA = ".";


        /// <summary>
        /// GTを作成する。
        /// </summary>
        /// <param name="gt">GT</param>
        public GT(string gt)
        {
            ArgumentException.ThrowIfNullOrEmpty(gt);

            Value = gt;
            IsNoData = gt == NO_DATA_GT;
            Indexes = IsNoData
                ? []
                : [.. gt.Split(_delimiter)
                    .Select(x => int.Parse(x))
                    .Distinct()];
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
        /// アレル塩基を取得する。
        /// </summary>
        /// <param name="allAlleles">全アレル塩基配列</param>
        /// <returns>アレル塩基</returns>
        public string ToAllele(string[] allAlleles)
        {
            if (allAlleles.Length == 0) throw new ArgumentException(null, nameof(allAlleles));

            if (IsNoData) return NO_DATA;

            var alleles = Indexes.Select(x => allAlleles[x]);
            return string.Join(ALLELE_DELIMITER, alleles);
        }
    }
}
